#include "camel/correct.h"

#include <cmath>
#include <fstream>
#include <functional>
#include <numeric>

#include "biosoup/timer.hpp"
#include "camel/coverage.h"
#include "camel/io.h"
#include "camel/mapping.h"
#include "detail/overlap.h"
#include "edlib.h"
#include "fmt/core.h"

namespace camel {

namespace detail {

static auto FindOverlapsAndFilterReads(
    State& state, std::vector<std::unique_ptr<biosoup::NucleicAcid>> src_reads)
    -> std::vector<ReadOverlapsPair> {
  auto reads_overlaps = std::vector<ReadOverlapsPair>();
  reads_overlaps.reserve(src_reads.size());

  auto overlaps = camel::FindOverlaps(state, MapCfg{}, src_reads);
  auto timer = biosoup::Timer();

  timer.Start();
  reads_overlaps.reserve(src_reads.size());
  std::transform(
      std::make_move_iterator(src_reads.begin()),
      std::make_move_iterator(src_reads.end()),
      std::make_move_iterator(overlaps.begin()),
      std::back_inserter(reads_overlaps),
      [](std::unique_ptr<biosoup::NucleicAcid> read,
         std::vector<biosoup::Overlap> ovlps) -> camel::ReadOverlapsPair {
        return {.read = std::move(read), .overlaps = std::move(ovlps)};
      });

  fmt::print(stderr,
             "[camel::detail::CollectOverlapsAndFilterReads]({:12.3f}) tied "
             "reads with overlaps\n",
             timer.Stop());

  {
    using namespace std::placeholders;
    using namespace std::literals;

    auto ovlp_cnsts = std::vector<std::uint64_t>(reads_overlaps.size());
    for (auto const& ro : reads_overlaps) {
      for (auto const& ovlp : ro.overlaps) {
        ++ovlp_cnsts[ovlp.lhs_id];
        ++ovlp_cnsts[ovlp.rhs_id];
      }
    }

    // preserving iterators
    auto unmapped_first = std::stable_partition(
        reads_overlaps.begin(), reads_overlaps.end(),
        [&ovlp_cnsts](camel::ReadOverlapsPair const& ro) -> bool {
          return ovlp_cnsts[ro.read->id] > 0UL;
        });

    auto const n_mapped = std::distance(reads_overlaps.begin(), unmapped_first);
    auto const n_unmapped = reads_overlaps.size() - n_mapped;

    timer.Start();

    // store unmapped reads
    if (unmapped_first != reads_overlaps.end()) {
      auto unmapped = std::vector<std::unique_ptr<biosoup::NucleicAcid>>();
      unmapped.reserve(n_unmapped);

      std::transform(std::make_move_iterator(unmapped_first),
                     std::make_move_iterator(reads_overlaps.end()),
                     std::back_inserter(unmapped),
                     [](camel::ReadOverlapsPair ro)
                         -> std::unique_ptr<biosoup::NucleicAcid> {
                       return std::move(ro.read);
                     });

      auto const dst_folder = state.log_path / "unmapped";
      camel::StoreSequences(state, unmapped, dst_folder);

      decltype(reads_overlaps)(std::make_move_iterator(reads_overlaps.begin()),
                               std::make_move_iterator(unmapped_first))
          .swap(reads_overlaps);
    }

    fmt::print(stderr,
               "[camel::CollectOverlapsAndFilterReads]({:12.3f}) stored {} / "
               "{} unmapped reads\n",
               timer.Stop(), n_unmapped, n_unmapped + n_mapped);
  }

  return reads_overlaps;
}

struct CmpOvlpLhsReadId {
  auto operator()(biosoup::Overlap const& ovlp,
                  std::uint32_t const read_id) const noexcept -> bool {
    return ovlp.lhs_id < read_id;
  }

  auto operator()(std::uint32_t const read_id,
                  biosoup::Overlap const& ovlp) const noexcept -> bool {
    return read_id < ovlp.lhs_id;
  }
};

static auto CollectOverlaps(State& state,
                            std::vector<ReadOverlapsPair> const& ros,
                            std::uint32_t const read_id)
    -> std::vector<biosoup::Overlap> {
  static constexpr auto kIntervalLen = 5000U;

  auto dst = std::vector<biosoup::Overlap>();
  auto const kNIntervals = std::ceil((1.0 * ros.size()) / kIntervalLen);

  auto interval_counts = std::vector<std::uint64_t>(kNIntervals);

  {
    auto count_futures = std::vector<std::future<std::uint64_t>>();

    count_futures.reserve(kNIntervals);
    for (auto first = ros.cbegin(); first != ros.cend();) {
      auto const last = std::min(std::next(first, kIntervalLen), ros.cend());

      count_futures.emplace_back(state.thread_pool->Submit(
          [read_id](std::vector<ReadOverlapsPair>::const_iterator first,
                    std::vector<ReadOverlapsPair>::const_iterator last)
              -> std::uint64_t {
            return std::transform_reduce(
                first, last, 0UL, std::plus<std::uint64_t>(),
                [read_id](ReadOverlapsPair const& ro) -> std::uint64_t {
                  auto const& ovlps = ro.overlaps;
                  if (!ovlps.empty() && ovlps.front().rhs_id == read_id) {
                    return ovlps.size();
                  } else {
                    auto const [lo, hi] =
                        std::equal_range(ovlps.cbegin(), ovlps.cend(), read_id,
                                         CmpOvlpLhsReadId());

                    return std::distance(lo, hi);
                  }
                });
          },
          first, last));

      first = last;
    }

    std::transform(
        count_futures.begin(), count_futures.end(), interval_counts.begin(),
        [](std::future<std::uint64_t>& f) -> std::uint64_t { return f.get(); });

    auto const n_total_ovlps =
        std::accumulate(interval_counts.cbegin(), interval_counts.cend(), 0UL);

    dst.resize(n_total_ovlps);
  }

  {
    auto collect_futures = std::vector<std::future<void>>();
    collect_futures.reserve(kNIntervals);

    auto interval_firsts = std::move(interval_counts);
    std::exclusive_scan(interval_firsts.begin(), interval_firsts.end(),
                        interval_firsts.begin(), 0UL);

    auto interval_id = 0U;
    for (auto first = ros.cbegin(); first != ros.cend(); ++interval_id) {
      auto const last = std::min(std::next(first, kIntervalLen), ros.cend());

      collect_futures.emplace_back(state.thread_pool->Submit(
          [read_id](std::vector<ReadOverlapsPair>::const_iterator first,
                    std::vector<ReadOverlapsPair>::const_iterator last,
                    std::vector<biosoup::Overlap>::iterator first_out) mutable
          -> void {
            for (; first != last; ++first) {
              auto const& ovlps = first->overlaps;
              if (!ovlps.empty() && ovlps.front().rhs_id == read_id) {
                std::transform(ovlps.begin(), ovlps.end(), first_out,
                               detail::ReverseOverlap);
                std::advance(first_out, ovlps.size());

              } else {
                auto const [lo, hi] = std::equal_range(
                    ovlps.cbegin(), ovlps.cend(), read_id, CmpOvlpLhsReadId());
                first_out = std::copy(lo, hi, first_out);
              }
            }
          },
          first, last, std::next(dst.begin(), interval_firsts[interval_id])));

      first = last;
    }

    std::for_each(collect_futures.begin(), collect_futures.end(),
                  std::mem_fn(&std::future<void>::wait));
  }

  return dst;
}

// TODO: extract to overlap detail
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

auto CalculateAlignments(State& state,
                         std::vector<ReadOverlapsPair> const& reads_overlaps,
                         std::vector<biosoup::Overlap> const& overlaps)
    -> std::vector<EdlibAlignResult> {
  auto dst = std::vector<EdlibAlignResult>();
  dst.reserve(overlaps.size());

  auto const calc_alignment =
      [&reads_overlaps](biosoup::Overlap const& ovlp) -> EdlibAlignResult {
    auto const [query_str, target_str] = OverlapStrings(reads_overlaps, ovlp);
    return edlibAlign(
        query_str.c_str(), query_str.size(), target_str.c_str(),
        target_str.size(),
        edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));
  };

  auto align_futures = std::vector<std::future<EdlibAlignResult>>();

  align_futures.reserve(overlaps.size());
  for (auto const& ovlp : overlaps) {
    align_futures.emplace_back(state.thread_pool->Submit(calc_alignment, ovlp));
  }

  std::transform(align_futures.begin(), align_futures.end(),
                 std::back_inserter(dst),
                 [](std::future<EdlibAlignResult>& f) -> EdlibAlignResult {
                   return f.get();
                 });

  return dst;
}

}  // namespace detail

auto SnpErrorCorrect(
    State& state, std::vector<std::unique_ptr<biosoup::NucleicAcid>> src_reads)
    -> std::vector<AnnotatedRead> {
  auto dst = std::vector<AnnotatedRead>();

  auto timer = biosoup::Timer();

  timer.Start();
  auto reads_overlaps =
      detail::FindOverlapsAndFilterReads(state, std::move(src_reads));
  fmt::print(stderr, "[camel]({:12.3f}) tied overlaps to reads\n",
             timer.Stop());

  timer.Start();
  auto const front_ovlps = detail::CollectOverlaps(
      state, reads_overlaps, reads_overlaps.back().read->id);
  fmt::print(stderr, "[camel]({:12.3f}) collected {} front overlaps\n",
             timer.Stop(), front_ovlps.size());

  timer.Start();
  auto const front_alignments =
      detail::CalculateAlignments(state, reads_overlaps, front_ovlps);
  fmt::print(stderr, "[camel]({:12.3f}) calculated front alignments\n", 
  timer.Stop());

  return dst;
}

}  // namespace camel
