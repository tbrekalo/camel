#include <concepts>
#include <cstdlib>
#include <fstream>
#include <future>
#include <optional>
#include <type_traits>

#include "camel/io.h"
#include "cxxopts.hpp"
#include "edlib.h"
#include "fmt/core.h"
#include "fmt/ranges.h"
#include "nlohmann/json.hpp"
#include "tbb/parallel_for.h"
#include "tbb/task_arena.h"
#include "tbb/task_group.h"
#include "tsl/robin_map.h"

using json = nlohmann::json;

static constexpr auto kMutPrefix = std::string_view("mut-");

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

struct ReadInfo {
  std::string origin_name;
  std::uint32_t origin_begin;
  std::uint32_t origin_length;
  bool is_muated;
  bool is_rc;
};

auto TokenViews(std::string_view src, std::string_view delim)
    -> std::vector<std::string_view> {
  auto dst = std::vector<std::string_view>();
  while (!src.empty()) {
    auto const delim_pos = src.find(delim);

    if (delim_pos > 0) {
      dst.push_back(src.substr(0, delim_pos));
      if (delim_pos > src.length()) {
        break;
      }
    }
    src.remove_prefix(delim_pos + delim.length());
  }

  return dst;
}

auto ParseUInt32(std::string_view src) {
  auto dst = std::uint32_t{0};
  for (auto digit : src) {
    dst = dst * 10 + (digit - '0');
  }
  return dst;
}

auto ParseReadInfo(std::string_view meta_data) -> ReadInfo {
  auto dst = ReadInfo{};
  if (meta_data.starts_with(kMutPrefix)) {
    dst.is_muated = true;

    meta_data = meta_data.substr(kMutPrefix.length());
    meta_data.remove_prefix(meta_data.find('-') + 1);
    meta_data.remove_prefix(meta_data.find('-') + 1);
  }

  dst.origin_name = meta_data.substr(0, meta_data.find('_'));
  meta_data.remove_prefix(meta_data.find('_'));

  auto tokens = TokenViews(meta_data, "_");
  dst.origin_begin = ParseUInt32(tokens[0]);
  dst.origin_length = ParseUInt32(tokens[5]);
  dst.is_rc = tokens[4] == "R";

  return dst;
};

struct EvalResult {
  std::uint64_t match_cnt;
  std::uint64_t mistmatch_cnt;
};

auto EvalRead(
    json const& mutations_index,
    tsl::robin_map<std::string, std::unique_ptr<biosoup::NucleicAcid>> const&
        references,
    std::unique_ptr<biosoup::NucleicAcid> read) -> EvalResult {
  auto dst = EvalResult{};

  auto read_info = ParseReadInfo(read->name);
  if (read_info.is_rc) {
    read->ReverseAndComplement();
  }

  auto const mutations =
      [](json const& m) -> std::vector<std::pair<std::uint32_t, char>> {
    auto dst = std::vector<std::pair<std::uint32_t, char>>();
    for (auto const& it : m.items()) {
      dst.emplace_back(ParseUInt32(it.key()), it.value().get<std::string>()[0]);
    }

    std::sort(dst.begin(), dst.end());
    return dst;
  }(mutations_index[read_info.origin_name]);

  auto target_str =
      references.at(read_info.origin_name)
          ->InflateData(read_info.origin_begin, read_info.origin_length);
  auto query_str = read->InflateData();

  auto edlib_res = edlibAlign(
      query_str.c_str(), query_str.size(), target_str.c_str(),
      target_str.size(),
      edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, nullptr, 0));

  auto target_pos = 0;
  auto query_pos = 0;

  auto j = std::distance(
      mutations.begin(),
      std::lower_bound(mutations.begin(), mutations.end(),
                       read_info.origin_begin,
                       [](auto const pic, auto pos) { return pic.first < pos; }

                       ));

  for (auto i = 0U; j < mutations.size() && i < edlib_res.alignmentLength;
       ++i) {
    auto const origin_pos = read_info.origin_begin + target_pos;
    if (mutations[j].first == origin_pos) {
      auto const is_match =
          query_str[query_pos] ==
          (read_info.is_muated ? mutations[j].second : target_str[target_pos]);

      dst.match_cnt += is_match;
      dst.mistmatch_cnt += !is_match;

      ++j;
    }

    query_pos += (edlib_res.alignment[i] != 2);
    target_pos += (edlib_res.alignment[i] != 1);
  }

  edlibFreeAlignResult(edlib_res);
  return dst;
}

int main(int argc, char** argv) {
  auto options =
      cxxopts::Options("eval", "tool for evaluating NanoSim corrected reads");

  /* clang-format off */
  options.add_options()
    ("q,query-reads", "reads to query against the reference",
      cxxopts::value<std::string>())
    ("r,reference", "reference genome",
      cxxopts::value<std::string>())
    ("j,json-mutations", "json file annoting mutations on the reference",
      cxxopts::value<std::string>())
    ("t,threads", "number of available threads for evaluation",
      cxxopts::value<std::uint32_t>()->default_value("1"))
    ("h,help", "print help");
  /* clang-format on */

  try {
    auto const result = options.parse(argc, argv);
    if (result.count("help")) {
      fmt::print(stderr, "{}", options.help());
      return EXIT_SUCCESS;
    }

    auto const n_threads = result["threads"].as<std::uint32_t>();
    auto const reads_path =
        std::filesystem::path(result["query-reads"].as<std::string>());
    auto const reference_path =
        std::filesystem::path(result["reference"].as<std::string>());
    auto const json_path = [&result]() -> std::optional<std::filesystem::path> {
      if (!result.count("json-mutations")) {
        return std::nullopt;
      }

      return std::filesystem::path(result["json-mutations"].as<std::string>());
    }();

    auto task_arena = tbb::task_arena(n_threads);
    task_arena.execute([&]() -> void {
      auto query_reads = std::vector<std::unique_ptr<biosoup::NucleicAcid>>();
      auto reference =
          tsl::robin_map<std::string, std::unique_ptr<biosoup::NucleicAcid>>();

      auto snps = std::optional<json>();

      {
        auto load_tg = tbb::task_group();
        auto const submit =
            [&load_tg]<class Fn, class... Args>(Fn && fn, Args && ... args)
          requires(std::invocable<Fn, Args...>)
        {
          using ReturnType = std::invoke_result_t<Fn, Args...>;
          auto task = std::make_unique<std::packaged_task<ReturnType()>>(
              [fn = std::forward<Fn>(fn),
               ... args = std::forward<Args>(args)]() mutable -> ReturnType {
                return std::forward<Fn>(fn)(std::forward<Args>(args)...);
              });

          auto future = task->get_future();
          load_tg.run([task = std::move(task)]() -> void { (*task)(); });

          return future;
        };

        auto query_reads_future =
            submit([&]() { return camel::LoadSequences(reads_path); });
        auto reference_future =
            submit([&]() { return camel::LoadSequences(reference_path); });

        auto snps_future = [&json_path,
                            &submit]() -> std::optional<std::future<json>> {
          if (!json_path) {
            return std::nullopt;
          }

          return submit(
              [](std::filesystem::path const& path) -> json {
                auto ifstrm = std::ifstream(path);
                return json::parse(ifstrm);
              },
              *json_path);
        }();

        query_reads = query_reads_future.get();
        for (auto&& it : reference_future.get()) {
          reference[it->name] = std::move(it);
        }

        if (snps_future) {
          snps = snps_future->get();
        }
      }

      auto eval_results = std::vector<EvalResult>(query_reads.size());
      tbb::parallel_for(
          std::size_t(0), query_reads.size(), [&](std::size_t i) -> void {
            eval_results[i] =
                EvalRead(*snps, reference, std::move(query_reads[i]));
          });

      auto matches = 0ULL;
      auto mismatches = 0ULL;

      for (auto const [m, n] : eval_results) {
        matches += m;
        mismatches += n;
      }

      fmt::print("{},{}", matches, mismatches);
    });

  } catch (std::exception const& e) {
    fmt::print(stderr, "{}", e.what());
  }

  return EXIT_SUCCESS;
}
