#include "camel/coverage.h"

#include <chrono>
#include <deque>
#include <functional>
#include <numeric>

#include "biosoup/timer.hpp"
#include "camel/io.h"
#include "detail/overlap.h"
#include "edlib.h"
#include "fmt/compile.h"
#include "fmt/core.h"

namespace camel {

namespace detail {

static constexpr std::size_t kAlignBatchCap = 1UL << 32UL;      // ~4gb
static constexpr std::size_t kDefaultSeqGroupSz = 1UL << 32UL;  // ~4.3gb

}  // namespace detail

auto CalculateCoverage(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool, MapCfg const map_cfg,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads,
    std::filesystem::path const& pile_storage_dir) -> void {
  auto read_overlap_vecs = FindOverlaps(thread_pool, map_cfg, reads);

  if (std::filesystem::exists(pile_storage_dir)) {
    std::filesystem::remove_all(pile_storage_dir);
  }

  auto find_batch_last =
      [&reads](
          std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator
              first,
          std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator
              last,
          std::size_t const batch_cap)
      -> std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator {
    for (auto batch_sz = 0UL; first != last && batch_sz < batch_cap; ++first) {
      batch_sz += (*first)->inflated_len * sizeof(Coverage);
    }

    return first;
  };

  std::filesystem::create_directory(pile_storage_dir);
  auto timer = biosoup::Timer();

  auto const overlap_strings =
      [&reads](
          biosoup::Overlap const& ovlp) -> std::pair<std::string, std::string> {
    auto query_str = reads[ovlp.lhs_id]->InflateData(
        ovlp.lhs_begin, ovlp.lhs_end - ovlp.lhs_begin);

    auto target_str = reads[ovlp.rhs_id]->InflateData(
        ovlp.rhs_begin, ovlp.rhs_end - ovlp.rhs_begin);

    if (!ovlp.strand) {
      auto rc = biosoup::NucleicAcid("", target_str);
      rc.ReverseAndComplement();

      target_str = rc.InflateData();
    }

    return std::pair(std::move(query_str), std::move(target_str));
  };

  auto const align_strings =
      [](std::string const& query_str,
         std::string const& target_str) -> EdlibAlignResult {
    return edlibAlign(
        query_str.c_str(), query_str.size(), target_str.c_str(),
        target_str.size(),
        edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));
  };

  auto const update_coverage = [&reads, &overlap_strings, &align_strings](
                                   Pile& pile,
                                   biosoup::Overlap const& ovlp) -> void {
    auto [query_substr, target_substr] = overlap_strings(ovlp);
    auto const edlib_res_align = align_strings(query_substr, target_substr);

    query_substr.clear();
    query_substr.shrink_to_fit();

    auto query_id = ovlp.lhs_id;
    auto target_id = ovlp.rhs_id;

    auto query_pos = ovlp.lhs_begin;
    auto target_pos = 0U;

    for (auto i = 0U; i < edlib_res_align.alignmentLength; ++i) {
      switch (edlib_res_align.alignment[i]) {
        case 0:
        case 3: {  // match
          /* clang-format off */
            switch (target_substr[target_pos]) {
              case 'A': ++pile.covgs[query_pos].a; break;
              case 'C': ++pile.covgs[query_pos].c; break;
              case 'G': ++pile.covgs[query_pos].g; break;
              case 'T': ++pile.covgs[query_pos].t; break;
              default: break;
            }
          /* clang-format on */

          ++query_pos;
          ++target_pos;

          break;
        }

        case 1: {  // insertion on target
          ++pile.covgs[query_pos].del;
          ++query_pos;
          break;
        }

        case 2: {  // insertion on query
          ++pile.covgs[query_pos].ins;
          ++target_id;
          break;
        }

        default: {
          break;
        }
      }
    }

    edlibFreeAlignResult(edlib_res_align);
  };

  auto coverage_futures = std::vector<std::future<Pile>>();

  auto const calc_coverage =
      [&read_overlap_vecs, &update_coverage](
          std::unique_ptr<biosoup::NucleicAcid> const& read) -> Pile {
    auto dst = Pile{.id = read->id,
                    .seq_name = read->name,
                    .covgs = std::vector<Coverage>(read->inflated_len)};

    auto const& read_ovlps = read_overlap_vecs[read->id];
    for (auto const& target_ovlp : read_ovlps) {
      auto const query_ovlp = detail::ReverseOverlap(target_ovlp);
      update_coverage(dst, query_ovlp);
    }

    for (auto remote_id = 0U; remote_id < read_overlap_vecs.size();
         ++remote_id) {
      if (remote_id != read->id) {
        auto const& remote_ovlps = read_overlap_vecs[remote_id];
        auto iter = std::lower_bound(
            remote_ovlps.cbegin(), remote_ovlps.cend(), read->id,
            [](biosoup::Overlap const& ovlp, std::uint32_t const query_id) {
              return ovlp.lhs_id < query_id;
            });

        for (; iter != remote_ovlps.cend() && iter->lhs_id == read->id;
             ++iter) {
          update_coverage(dst, *iter);
        }
      }
    }

    return dst;
  };

  auto const calc_coverage_async =
      [&thread_pool, &calc_coverage](
          std::reference_wrapper<std::unique_ptr<biosoup::NucleicAcid> const>
              read) -> std::future<Pile> {
    return thread_pool->Submit(calc_coverage, read);
  };

  auto serialize_futures = std::deque<std::future<
      std::tuple<std::filesystem::path, std::size_t, std::size_t, double>>>();

  auto is_ser_future_ready =
      [](decltype(serialize_futures)::const_reference f) -> bool {
    return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
  };

  auto pop_and_report = [&serialize_futures]() -> void {
    auto const [dst_file, raw_sz, com_sz, comp_time] =
        serialize_futures.front().get();
    serialize_futures.pop_front();

    fmt::print(stderr,
               FMT_COMPILE("[camel::CalculateCoverage]({:12.3f}) {} -> "
                           "comp ratio: {:1.3f}\n"),
               comp_time, dst_file.string(), 1. * raw_sz / com_sz);
  };

  {
    auto batch_id = 0U;
    auto active_piles = std::vector<Pile>();
    for (auto align_first = reads.cbegin(); align_first != reads.cend();
         ++batch_id) {
      timer.Start();
      auto const align_last =
          find_batch_last(align_first, reads.cend(), detail::kAlignBatchCap);

      auto const batch_sz = std::distance(align_first, align_last);
      coverage_futures.reserve(batch_sz);
      active_piles.reserve(batch_sz);

      while (!serialize_futures.empty() &&
             is_ser_future_ready(serialize_futures.front())) {
        pop_and_report();
      }

      std::transform(align_first, align_last,
                     std::back_inserter(coverage_futures), calc_coverage_async);

      fmt::print(
          stderr,
          "[camel::CalculateCoverage] calculating coverage for {} reads\n",
          batch_sz);

      std::transform(coverage_futures.begin(), coverage_futures.end(),
                     std::back_inserter(active_piles),
                     std::mem_fn(&std::future<Pile>::get));

      serialize_futures.emplace_back(thread_pool->Submit(
          [&pile_storage_dir](std::vector<Pile> batch_piles,
                              std::size_t batch_id)
              -> std::tuple<std::filesystem::path, std::size_t, std::size_t,
                            double> {
            auto ser_timer = biosoup::Timer();

            ser_timer.Start();
            auto const dst_file = SerializePileBatch(
                batch_piles.cbegin(), batch_piles.cend(), pile_storage_dir,
                fmt::format("pile_batch_{:04d}", batch_id));

            auto const uncompressed_bytes = std::transform_reduce(
                std::make_move_iterator(batch_piles.begin()),
                std::make_move_iterator(batch_piles.end()), 0UL,
                std::plus<std::size_t>(), [](Pile&& pile) -> std::size_t {
                  return sizeof(Pile::id) + pile.seq_name.size() +
                         pile.covgs.size() * sizeof(Coverage);
                });

            auto const compressed_bytes = std::filesystem::file_size(dst_file);

            return std::make_tuple(std::move(dst_file), uncompressed_bytes,
                                   compressed_bytes, ser_timer.Stop());
          },
          std::move(active_piles), batch_id));

      fmt::print(
          stderr,
          "[camel::CalculateCoverage]({:12.3f}) finished coverage batch\n",
          timer.Stop(), batch_sz);

      coverage_futures.clear();
      align_first = align_last;
    }
  }

  timer.Start();
  while (!serialize_futures.empty()) {
    pop_and_report();
  }
  timer.Stop();

  fmt::print(stderr, "[camel::CalculateCoverage]({:12.3f})\n",
             timer.elapsed_time());
}

}  // namespace camel
