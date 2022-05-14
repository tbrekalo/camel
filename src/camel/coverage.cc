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
#include "fmt/ostream.h"

namespace camel {

namespace detail {

static constexpr std::size_t kAlignBatchCap = 1UL << 32UL;      // ~4gb
static constexpr std::size_t kDefaultSeqGroupSz = 1UL << 32UL;  // ~4.3gb

struct InsSignal {
  using ValueType = std::uint16_t;

  ValueType a;
  ValueType c;
  ValueType t;
  ValueType g;
};

[[nodiscard]] static auto CmpReadsOvlpsPairsByReadId(
    ReadOverlapsPair const& lhs, ReadOverlapsPair const& rhs) -> bool {
  return lhs.read->id < rhs.read->id;
}

[[nodiscard]] static auto FindBatchLast(
    std::vector<ReadOverlapsPair>::const_iterator first,
    std::vector<ReadOverlapsPair>::const_iterator last, std::size_t batch_cap) {
  for (auto batch_sz = 0UL; first != last && batch_sz < batch_cap; ++first) {
    batch_sz += first->read->inflated_len * sizeof(Coverage);
  }

  return first;
}

[[nodiscard]] static auto OverlapStrings(
    std::vector<ReadOverlapsPair> const& reads_overlaps,
    biosoup::Overlap const& ovlp) -> std::pair<std::string, std::string> {
  auto const query_read_iter = std::lower_bound(
      reads_overlaps.cbegin(), reads_overlaps.cend(), ovlp.lhs_id,
      [](ReadOverlapsPair const& ro, std::uint32_t const query_id) -> bool {
        return ro.read->id < query_id;
      });

  auto const target_read_iter = std::lower_bound(
      reads_overlaps.cbegin(), reads_overlaps.cend(), ovlp.rhs_id,
      [](ReadOverlapsPair const& ro, std::uint32_t const target_id) -> bool {
        return ro.read->id < target_id;
      });

  auto query_str = query_read_iter->read->InflateData(
      ovlp.lhs_begin, ovlp.lhs_end - ovlp.lhs_begin);

  auto target_str = target_read_iter->read->InflateData(
      ovlp.rhs_begin, ovlp.rhs_end - ovlp.rhs_begin);

  if (!ovlp.strand) {
    auto rc = biosoup::NucleicAcid("", target_str);
    rc.ReverseAndComplement();

    target_str = rc.InflateData();
  }

  return std::pair(std::move(query_str), std::move(target_str));
}

[[nodiscard]] static auto AlignStrings(std::string const& query_str,
                                       std::string const& target_str)
    -> EdlibAlignResult {
  return edlibAlign(
      query_str.c_str(), query_str.size(), target_str.c_str(),
      target_str.size(),
      edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));
}

}  // namespace detail

CAMEL_EXPORT auto CalculateCoverage(
    State& state, std::vector<ReadOverlapsPair> const& reads_overlaps,
    std::filesystem::path const& pile_storage_dir) -> void {
  // consider making reads_overlaps just &; not const&
  if (!std::is_sorted(reads_overlaps.cbegin(), reads_overlaps.cend(),
                      detail::CmpReadsOvlpsPairsByReadId)) {
    throw std::runtime_error(
        "[camel::CalculateCoverage] input is exptected to be sorted");
  }

  if (std::filesystem::exists(pile_storage_dir)) {
    std::filesystem::remove_all(pile_storage_dir);
  }

  std::filesystem::create_directory(pile_storage_dir);
  auto timer = biosoup::Timer();

  auto const update_coverage =
      [&reads_overlaps](Pile& pile, biosoup::Overlap const& ovlp) -> void {
    auto [query_substr, target_substr] =
        detail::OverlapStrings(reads_overlaps, ovlp);
    auto const edlib_res_align =
        detail::AlignStrings(query_substr, target_substr);

    query_substr.clear();
    query_substr.shrink_to_fit();

    auto query_pos = ovlp.lhs_begin;
    auto target_pos = 0U;

    for (auto i = 0U; i < edlib_res_align.alignmentLength; ++i) {
      switch (edlib_res_align.alignment[i]) {
        case 0:
        case 3: {  // mismatch
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
          ++target_pos;
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
      [&reads_overlaps,
       &update_coverage](ReadOverlapsPair const& read_overlaps) -> Pile {
    auto dst =
        Pile{.id = read_overlaps.read->id,
             .seq_name = read_overlaps.read->name,
             .covgs = std::vector<Coverage>(read_overlaps.read->inflated_len)};

    auto const& local_ovlps = read_overlaps.overlaps;
    for (auto const& target_ovlp : local_ovlps) {
      auto const query_ovlp = detail::ReverseOverlap(target_ovlp);
      update_coverage(dst, query_ovlp);
    }

    for (auto const& [other_read, remote_ovlps] : reads_overlaps) {
      if (read_overlaps.read->id != other_read->id) {
        auto iter = std::lower_bound(
            remote_ovlps.cbegin(), remote_ovlps.cend(), read_overlaps.read->id,
            [](biosoup::Overlap const& ovlp, std::uint32_t const query_id) {
              return ovlp.lhs_id < query_id;
            });

        for (; iter != remote_ovlps.cend() &&
               iter->lhs_id == read_overlaps.read->id;
             ++iter) {
          update_coverage(dst, *iter);
        }
      }
    }

    return dst;
  };

  auto const calc_coverage_async =
      [&thread_pool = state.thread_pool, &calc_coverage](
          std::reference_wrapper<ReadOverlapsPair const> read_overlaps)
      -> std::future<Pile> {
    return thread_pool->Submit(calc_coverage, read_overlaps);
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
    for (auto align_first = reads_overlaps.cbegin();
         align_first != reads_overlaps.cend(); ++batch_id) {
      timer.Start();
      auto const align_last = detail::FindBatchLast(
          align_first, reads_overlaps.cend(), detail::kAlignBatchCap);

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

      serialize_futures.emplace_back(state.thread_pool->Submit(
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

      fmt::print(stderr,
                 "[camel::CalculateCoverage]({:12.3f}) finished coverage "
                 "batch; {:3.2f}\% reads covered\n",
                 timer.Stop(),
                 100. * std::distance(reads_overlaps.cbegin(), align_last) /
                     reads_overlaps.size());

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
