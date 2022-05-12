#include <algorithm>
#include <array>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

#include "biosoup/timer.hpp"
#include "camel/correct.h"
#include "camel/coverage.h"
#include "camel/detail/overlap.h"
#include "camel/io.h"
#include "camel/mapping.h"
#include "camel/state.h"
#include "cxxopts.hpp"
#include "fmt/core.h"
#include "fmt/ostream.h"
#include "jemalloc/jemalloc.h"
#include "thread_pool/thread_pool.hpp"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects(0);

auto main(int argc, char** argv) -> int {
  auto options =
      cxxopts::Options("camel", "Camel is haplotype aware detection tool");

  /* clang-format off */
  options.add_options()
    ("t,threads", "number of threads avalable for execution",
            cxxopts::value<std::uint32_t>())
    ("s,serialization_dst", "destination folder for pile serialization",
            cxxopts::value<std::string>()->default_value("./camel_piles"))
    ("l,log_dst", "destination for for logging information",
            cxxopts::value<std::string>()->default_value("./camel_log"))
    ("paths", "input fastq reads", 
            cxxopts::value<std::vector<std::string>>());
  options.parse_positional("paths");
  /* clang-format on */

  auto const result = options.parse(argc, argv);

  auto state =
      camel::State{.thread_pool = std::make_shared<thread_pool::ThreadPool>(
                       result["threads"].as<std::uint32_t>()),
                   .log_path = result["log_dst"].as<std::string>()};

  auto const ser_dst_path =
      std::filesystem::path(result["serialization_dst"].as<std::string>());

  if (std::filesystem::exists(state.log_path)) {
    std::filesystem::remove_all(state.log_path);
  }

  std::filesystem::create_directory(state.log_path);

  auto paths = std::vector<std::filesystem::path>();

  {
    auto paths_strs = result["paths"].as<std::vector<std::string>>();
    for (auto const& it : paths_strs) {
      auto it_path = std::filesystem::path(std::move(it));
      if (std::filesystem::is_regular_file(it_path)) {
        paths.emplace_back(std::move(it_path));
      } else if (std::filesystem::is_directory(it_path)) {
        for (auto dir_entry : std::filesystem::recursive_directory_iterator(
                 std::move(it_path))) {
          if (std::filesystem::is_regular_file(dir_entry)) {
            paths.emplace_back(std::move(dir_entry));
          }
        }
      }
    }
  }

  decltype(paths)(std::make_move_iterator(paths.begin()),
                  std::make_move_iterator(paths.end()))
      .swap(paths);

  auto timer = biosoup::Timer();
  auto reads_overlaps = std::vector<camel::ReadOverlapsPair>();

  {
    timer.Start();
    auto reads = camel::LoadSequences(state, paths);
    fmt::print(stderr, "[camel]({:12.3f}) loaded {} reads\n", timer.Stop(),
               reads.size());

    timer.Start();
    auto overlaps = camel::FindConfidentOverlaps(state, camel::MapCfg{}, reads);

    auto ovlp_futures = std::vector<std::future<void>>();
    ovlp_futures.reserve(reads.size());

    auto const transform_target =
        [&overlaps](std::uint32_t const target_id) -> void {
      auto const target_cnt = std::count_if(
          overlaps[target_id].cbegin(), overlaps[target_id].cend(),
          [target_id](biosoup::Overlap const& ovlp) -> bool {
            return ovlp.rhs_id < target_id;
          });

      auto dst = std::vector<biosoup::Overlap>();
      dst.reserve(target_cnt);

      for (auto const& ovlp : overlaps[target_id]) {
        if (ovlp.rhs_id < target_id) {
          dst.push_back(camel::detail::ReverseOverlap(ovlp));
        }
      }

      auto const by_query_id_cmp = [](biosoup::Overlap const& a,
                                      biosoup::Overlap const& b) -> bool {
        return a.lhs_id < b.lhs_id;
      };

      std::sort(dst.begin(), dst.end(), by_query_id_cmp);
      std::swap(overlaps[target_id], dst);
    };

    auto const transform_range = [&transform_target](
                                     std::uint32_t const first_id,
                                     std::uint32_t const last_id) -> void {
      for (auto target_id = first_id; target_id != last_id; ++target_id) {
        transform_target(target_id);
      }
    };

    for (auto target_id = 0U; target_id < reads.size(); target_id += 1000U) {
      ovlp_futures.emplace_back(state.thread_pool->Submit(
          transform_range, target_id,
          std::min(target_id + 1000U,
                   static_cast<std::uint32_t>(reads.size()))));
    }

    // for (auto target_id = 0U; target_id < reads.size(); ++target_id) {
    //   ovlp_futures.emplace_back(
    //       thread_pool->Submit(transform_target, target_id));
    // }

    std::for_each(ovlp_futures.begin(), ovlp_futures.end(),
                  std::mem_fn(&std::future<void>::get));

    // auto const n_ovlps = std::transform_reduce(
    //     overlaps.cbegin(), overlaps.cend(), 0UL, std::plus<std::uint64_t>(),
    //     [](std::vector<biosoup::Overlap> const& ovlp_vec) -> std::uint64_t {
    //       return ovlp_vec.size();
    //     });

    // fmt::print(stderr, "[camel::de] {}\n", n_ovlps);

    reads_overlaps.reserve(reads.size());
    std::transform(
        std::make_move_iterator(reads.begin()),
        std::make_move_iterator(reads.end()),
        std::make_move_iterator(overlaps.begin()),
        std::back_inserter(reads_overlaps),
        [](std::unique_ptr<biosoup::NucleicAcid> read,
           std::vector<biosoup::Overlap> ovlps) -> camel::ReadOverlapsPair {
          return {.read = std::move(read), .overlaps = std::move(ovlps)};
        });

    fmt::print(stderr, "[camel]({:12.3f}) tied reads with overlaps\n",
               timer.Stop());
  }

  {
    timer.Start();

    using namespace std::placeholders;
    using namespace std::literals;

    auto ovlp_cnsts = std::vector<std::uint64_t>(reads_overlaps.size());
    for (auto const& ro : reads_overlaps) {
      for (auto const& ovlp : ro.overlaps) {
        ++ovlp_cnsts[ovlp.lhs_id];
        ++ovlp_cnsts[ovlp.rhs_id];
      }
    }

    auto unmapped_first = std::stable_partition(
        reads_overlaps.begin(), reads_overlaps.end(),
        [&ovlp_cnsts](camel::ReadOverlapsPair const& ro) -> bool {
          return ovlp_cnsts[ro.read->id] > 0UL;
        });

    auto const n_mapped = std::distance(reads_overlaps.begin(), unmapped_first);
    auto const n_unmapped = reads_overlaps.size() - n_mapped;

    // store unmapped reads
    if (unmapped_first != reads_overlaps.end()) {
      auto const is_fasta = unmapped_first->read->block_quality.empty();

      auto const unmapped_files_log =
          state.log_path / ("unmapped"s + (is_fasta ? ".fa"s : ".fq"s));

      auto ofstrm =
          std::fstream(unmapped_files_log, std::ios::out | std::ios::trunc);

      auto const store_fasta_impl =
          +[](std::ostream& ostrm, camel::ReadOverlapsPair const& ro) -> void {
        ostrm << '>' << ro.read->name << '\n' << ro.read->InflateData() << '\n';
      };

      auto const store_fastq_impl =
          +[](std::ostream& ostrm, camel::ReadOverlapsPair const& ro) -> void {
        ostrm << '@' << ro.read->name << '\n'
              << ro.read->InflateData() << '\n'
              << "+\n"
              << ro.read->InflateQuality() << '\n';
      };

      using StoreFnSig = void(std::ostream&, camel::ReadOverlapsPair const&);
      using StoreFnPtr = std::add_const_t<std::add_pointer_t<StoreFnSig>>;

      static_assert(std::is_same_v<decltype(store_fasta_impl), StoreFnPtr>);
      static_assert(std::is_same_v<decltype(store_fastq_impl), StoreFnPtr>);

      auto const store_fn_impl = is_fasta ? store_fasta_impl : store_fastq_impl;

      auto const store_fn = [&ofstrm, impl = store_fn_impl](
                                camel::ReadOverlapsPair const& ro) -> void {
        impl(ofstrm, ro);
      };

      std::for_each(std::make_move_iterator(unmapped_first),
                    std::make_move_iterator(reads_overlaps.end()), store_fn);

      reads_overlaps.erase(unmapped_first, reads_overlaps.end());
      decltype(reads_overlaps)(std::make_move_iterator(reads_overlaps.begin()),
                               std::make_move_iterator(reads_overlaps.end()))
          .swap(reads_overlaps);

      fmt::print(stderr, "[camel]({:12.3f}) stored {} unmapped reads\n",
                 timer.Stop(), n_unmapped);
    }

    {
      timer.Start();
      camel::CalculateCoverage(state, reads_overlaps, ser_dst_path);
      timer.Stop();
    }
  }

  fmt::print(stderr, "[camel]({:12.3f}) done\n", timer.elapsed_time());

  return EXIT_SUCCESS;
}
