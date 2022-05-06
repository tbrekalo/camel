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
#include "cxxopts.hpp"
#include "fmt/core.h"
#include "fmt/ostream.h"
// #include "jemalloc/jemalloc.h"
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

  auto const thread_pool = std::make_shared<thread_pool::ThreadPool>(
      result["threads"].as<std::uint32_t>());

  auto const ser_dst_path =
      std::filesystem::path(result["serialization_dst"].as<std::string>());

  auto const log_dst_path =
      std::filesystem::path(result["log_dst"].as<std::string>());

  if (std::filesystem::exists(log_dst_path)) {
    std::filesystem::remove_all(log_dst_path);
  }

  std::filesystem::create_directory(log_dst_path);

  auto paths = std::vector<std::filesystem::path>();
  auto paths_strs = result["paths"].as<std::vector<std::string>>();
  std::transform(std::make_move_iterator(paths_strs.begin()),
                 std::make_move_iterator(paths_strs.end()),
                 std::back_inserter(paths),
                 [](std::string&& path_str) -> std::filesystem::path {
                   return std::filesystem::path(std::move(path_str));
                 });

  auto timer = biosoup::Timer();

  auto reads_overlaps = std::vector<camel::ReadOverlapsPair>();

  // {
  //   timer.Start();

  //   constexpr auto kExpectedFileSize = 1UL << 26UL;  // 64MiB
  //   auto reads = camel::LoadSequences(thread_pool, paths);
  //   fmt::print(stderr, "[camel]({:12.3f}) loaded {} reads\n", timer.Stop(),
  //              reads.size());

  //   timer.Start();

  //   auto const dst_dir_path =
  //       std::filesystem::path("/storage2/tbrekalo/yeast_reads/");
  //   auto dst_futures = std::vector<std::future<void>>();

  //   auto const find_batch_last =
  //       [kExpectedFileSize](
  //           decltype(reads)::iterator first,
  //           decltype(reads)::iterator last) -> decltype(reads)::iterator {
  //     for (auto batch_sz = 0UL; batch_sz < kExpectedFileSize && first !=
  //     last;
  //          ++first) {
  //       batch_sz += (*first)->inflated_len * 2 + (*first)->name.size();
  //     }

  //     return first;
  //   };

  //   auto batch_id = 0U;
  //   for (auto batch_first = reads.begin(); batch_first != reads.end();
  //        ++batch_id) {
  //     auto const batch_last = find_batch_last(batch_first, reads.end());

  //     fmt::print(stderr, "[de] batch_{:04d} - n seqs -> {}\n", batch_id,
  //                std::distance(batch_first, batch_last));

  //     dst_futures.emplace_back(thread_pool->Submit(
  //         [&dst_dir_path](auto first, auto last,
  //                         std::uint32_t const batch_id) -> void {
  //           auto const dst_path =
  //               dst_dir_path /
  //               fmt::format("saccharomyces_cerevisiae_{:04d}.fastq",
  //               batch_id);
  //           auto dst_strm =
  //               std::fstream(dst_path, std::ios::out | std::ios::trunc);

  //           for (; first != last; ++first) {
  //             fmt::print(dst_strm, "@{}\n{}\n+\n{}\n", (*first)->name,
  //                        (*first)->InflateData(),
  //                        (*first)->InflateQuality());
  //           }
  //         },
  //         std::make_move_iterator(batch_first),
  //         std::make_move_iterator(batch_last), batch_id));

  //     batch_first = batch_last;
  //   }

  //   std::for_each(dst_futures.begin(), dst_futures.end(),
  //                 std::mem_fn(&std::future<void>::wait));

  //   fmt::print(stderr, "[camel]({:12.3f}) distributed reads\n",
  //   timer.Stop());
  // }

  {
    timer.Start();
    auto reads = camel::LoadSequences(thread_pool, paths);
    fmt::print(stderr, "[camel]({:12.3f}) loaded {} reads\n", timer.Stop(),
               reads.size());

    // {
    //   timer.Start();
    //   auto res = camel::SnpErrorCorrect(thread_pool, std::move(reads));
    //   timer.Stop();
    // }
    timer.Start();
    auto overlaps =
        camel::FindConfidentOverlaps(thread_pool, camel::MapCfg{}, reads);

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

    for (auto target_id = 0U; target_id < reads.size(); ++target_id) {
      ovlp_futures.emplace_back(
          thread_pool->Submit(transform_target, target_id));
    }

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
          log_dst_path / ("unmapped"s + (is_fasta ? ".fa"s : ".fq"s));

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
      camel::CalculateCoverage(thread_pool, reads_overlaps, ser_dst_path);
      timer.Stop();
    }
  }

  fmt::print(stderr, "[camel]({:12.3f}) done\n", timer.elapsed_time());

  return EXIT_SUCCESS;
}
