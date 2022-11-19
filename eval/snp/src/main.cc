#include <concepts>
#include <cstdlib>
#include <fstream>
#include <future>
#include <optional>
#include <type_traits>

#include "camel/io.h"
#include "cxxopts.hpp"
#include "fmt/core.h"
#include "nlohmann/json.hpp"
#include "tbb/task_arena.h"
#include "tbb/task_group.h"

using json = nlohmann::json;

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

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
      auto reference = std::vector<std::unique_ptr<biosoup::NucleicAcid>>();
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
        reference = reference_future.get();
        if (snps_future) {
          snps = snps_future->get();
        }

        fmt::print(stderr, "[snp_eval] loaded {} query reads\n",
                   query_reads.size());
        fmt::print(stderr, "[snp_eval] loaded {} reference contigs\n",
                   reference.size());
      }
    });

  } catch (std::exception const& e) {
    fmt::print(stderr, "{}", e.what());
  }

  return EXIT_SUCCESS;
}
