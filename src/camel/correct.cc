#include "camel/correct.h"

#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <numeric>
#include <random>
#include <variant>

#include "biosoup/timer.hpp"
#include "camel/io.h"
#include "camel/mapping.h"
#include "detail/overlap.h"
#include "edlib.h"
#include "fmt/core.h"
#include "fmt/ostream.h"
#include "fmt/ranges.h"
#include "ram/minimizer_engine.hpp"
#include "spoa/alignment_engine.hpp"
#include "spoa/spoa.hpp"
#include "tsl/robin_map.h"

namespace camel {

namespace detail {

static auto constexpr kSqrt5 = 2.2360679775;
static auto constexpr kShrinkShift = 3U;
static auto constexpr kWinLength = 420U;

struct CoverageSignals {
  std::array<std::uint16_t, 6U> signals;

  static constexpr auto kDelIdx = 4U;
  static constexpr auto kInsIdx = 5U;
};

struct Interval {
  std::uint32_t first;
  std::uint32_t last;
};

constexpr auto operator==(Interval const lhs, Interval const rhs) noexcept
    -> bool {
  return lhs.first == rhs.first && lhs.last == rhs.last;
}

constexpr auto operator!=(Interval const lhs, Interval const rhs) -> bool {
  return !(lhs == rhs);
}

constexpr auto IntervalLength(Interval const intv) -> std::uint32_t {
  return intv.last - intv.first;
}

struct AlignedSegment {
  Interval aligned_interval;
  std::string bases;
};

struct ReferenceWindow {
  Interval interval;
  std::vector<AlignedSegment> aligned_segments;
};

class NucleicView {
 public:
  NucleicView(biosoup::NucleicAcid const* nucleic_acid,
              bool is_reverse_complement)
      : nucleic_acid_(nucleic_acid),
        fetch_code_impl_(is_reverse_complement ? &FetchReverseComplementCodeImpl
                                               : &FetchCodeImpl) {}

  auto Code(std::size_t const pos) const noexcept -> std::uint8_t {
    return fetch_code_impl_(nucleic_acid_, pos);
  }

  auto InflateData(std::uint32_t const pos, std::uint32_t const len) const
      -> std::string {
    auto dst = std::string(len, '\0');
    for (auto i = 0U; i < len; ++i) {
      dst[i] = biosoup::kNucleotideDecoder[Code(pos + i)];
    }
    return dst;
  }

 private:
  using FetchCodeImplPtr = std::uint8_t (*)(biosoup::NucleicAcid const*,
                                            std::size_t);

  static auto FetchCodeImpl(biosoup::NucleicAcid const* nucleic_acid,
                            std::size_t pos) noexcept -> std::uint8_t {
    return ((nucleic_acid->deflated_data[pos >> 5] >> ((pos << 1) & 63)) & 3);
  }

  static auto FetchReverseComplementCodeImpl(
      biosoup::NucleicAcid const* nucleic_acid, std::size_t pos) noexcept
      -> std::uint8_t {
    return FetchCodeImpl(nucleic_acid, nucleic_acid->inflated_len - 1 - pos) ^
           3;
  }

  biosoup::NucleicAcid const* nucleic_acid_;
  FetchCodeImplPtr fetch_code_impl_;
};

[[nodiscard]] static auto EstimateCoverage(
    State& state,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads,
    std::vector<std::vector<biosoup::Overlap>> const& overlaps)
    -> std::uint32_t {
  auto covg_futures = std::vector<std::future<std::uint16_t>>();
  covg_futures.reserve(reads.size());

  auto const covg_task =
      [&reads, &overlaps](std::uint32_t const read_id) -> std::uint16_t {
    auto pile = std::vector<std::uint16_t>(
        (reads[read_id]->inflated_len >> kShrinkShift) + 2U);

    auto events = std::vector<std::uint32_t>();
    events.reserve(2 * overlaps[read_id].size());

    for (auto const& ovlp : overlaps[read_id]) {
      events.push_back(((ovlp.lhs_begin >> kShrinkShift) + 1U) << 1U);
      events.push_back((((ovlp.lhs_end >> kShrinkShift) - 1U) << 1U) | 1U);
    }

    auto covg = 0U;
    auto last_event = 0U;
    std::sort(events.begin(), events.end());
    for (auto const event : events) {
      for (auto i = last_event; i < (event >> 1U); ++i) {
        pile[i] = std::clamp(pile[i] + covg, 0U,
                             1U * std::numeric_limits<std::uint16_t>::max());
      }

      covg += 1 - (2 * (event & 1));
      last_event = event >> 1U;
    }

    std::nth_element(pile.begin(), pile.begin() + pile.size() / 2, pile.end());
    return pile[pile.size() / 2];
  };

  for (auto read_id = 0U; read_id < reads.size(); ++read_id) {
    if (!overlaps[read_id].empty()) {
      covg_futures.emplace_back(state.thread_pool->Submit(covg_task, read_id));
    }
  }

  auto covgs = std::vector<std::uint16_t>(covg_futures.size());
  std::transform(covg_futures.begin(), covg_futures.end(), covgs.begin(),
                 std::mem_fn(&std::future<std::uint16_t>::get));

  std::nth_element(covgs.begin(), covgs.begin() + covgs.size() / 2,
                   covgs.end());

  return covgs[covgs.size() / 2];
}

[[nodiscard]] static auto ExtractSubstrings(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads,
    biosoup::Overlap const& ovlp) -> std::pair<std::string, std::string> {
  auto const lhs_view = NucleicView(reads[ovlp.lhs_id].get(), false);
  auto const rhs_view = NucleicView(reads[ovlp.rhs_id].get(),
                                    /*is_reverse_complement = */ !ovlp.strand);

  auto lhs_str =
      lhs_view.InflateData(ovlp.lhs_begin, ovlp.lhs_end - ovlp.lhs_begin);
  auto rhs_str =
      rhs_view.InflateData(ovlp.rhs_begin, ovlp.rhs_end - ovlp.rhs_begin);

  return {std::move(lhs_str), std::move(rhs_str)};
}

[[nodiscard]] static auto AlignStrings(std::string_view lhs_str_view,
                                       std::string_view rhs_str_view)
    -> EdlibAlignResult {
  /* clang-format off */
  return edlibAlign(
      lhs_str_view.data(), lhs_str_view.length(), 
      rhs_str_view.data(), rhs_str_view.length(),
      edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));
  /* clang-format on */
}

[[nodiscard]] static auto CalculateCoverage(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads,
    std::vector<biosoup::Overlap> const& overlaps,
    std::vector<EdlibAlignResult> const& edlib_results)
    -> std::vector<CoverageSignals> {
  auto const query_id = overlaps.front().lhs_id;
  auto dst = std::vector<CoverageSignals>(reads[query_id]->inflated_len + 1U);

  for (auto i = 0U; i < reads[query_id]->inflated_len; ++i) {
    ++dst[i].signals[reads[query_id]->Code(i)];
  }

  for (auto ovlp_idx = 0U; ovlp_idx < overlaps.size(); ++ovlp_idx) {
    auto const& edlib_res = edlib_results[ovlp_idx];
    auto const rhs_view = NucleicView(reads[overlaps[ovlp_idx].rhs_id].get(),
                                      !overlaps[ovlp_idx].strand);

    auto lhs_pos = overlaps[ovlp_idx].lhs_begin;
    auto rhs_pos = overlaps[ovlp_idx].rhs_begin;
    for (auto i = 0; i < edlib_res.alignmentLength; ++i) {
      // if (lhs_pos >= reads[query_id]->inflated_len) {
      //   fmt::print(stderr, "({} >= {}) ... (alignment[{} / {}] = {})\n",
      //   lhs_pos,
      //              reads[query_id]->inflated_len, i,
      //              edlib_res.alignmentLength,
      //              static_cast<int>(edlib_res.alignment[i]));
      // }
      switch (edlib_res.alignment[i]) {
        case 0:  // match
        case 3:  // mismatch
          ++dst[lhs_pos].signals[rhs_view.Code(rhs_pos)];
          ++lhs_pos;
          ++rhs_pos;
          break;
        case 1:  // deletion on the query
          ++dst[lhs_pos].signals[CoverageSignals::kDelIdx];
          ++lhs_pos;
          break;
        case 2:  // insertion on the query
          // ++dst[lhs_pos].signals[CoverageSignals::kInsIdx];
          ++rhs_pos;
          break;
        default:
          break;
      }
    }
  }

  return dst;
}

[[nodiscard]] static auto IsUnstableSite(CoverageSignals const& covg,
                                         std::uint32_t const covg_estimate,
                                         double const strong_base_ratio,
                                         double const indle_ratio) -> bool {
  auto const sig_sum =
      std::accumulate(covg.signals.cbegin(), covg.signals.cend(), 0UL,
                      std::plus<std::uint32_t>());
  auto const strong_base_threshold =
      static_cast<std::uint32_t>(std::round(strong_base_ratio * sig_sum));
  auto const indle_threshold =
      static_cast<std::uint32_t>(std::round(indle_ratio * sig_sum));
  auto const cutoff_threshold =
      static_cast<std::uint32_t>(std::round(covg_estimate * kSqrt5));

  auto dst = true;
  for (auto i = 0U; dst && i < 4U; ++i) {
    if (covg.signals[i] >= strong_base_threshold ||
        covg.signals[i] < cutoff_threshold) {
      dst = false;
    }
  }

  for (auto i = 4U; dst && i < 6U; ++i) {
    if (covg.signals[i] > indle_ratio) {
      dst = false;
    }
  }

  return dst;
}

[[nodiscard]] static auto IsSnpSite(CoverageSignals const& covg,
                                    std::uint32_t const covg_estimate) -> bool {
  return IsUnstableSite(covg, covg_estimate, 0.8, 0.4);
}

[[nodiscard]] static auto IsErrorSite(CoverageSignals const& covg,
                                      std::uint32_t const covg_estimate)
    -> bool {
  return IsUnstableSite(covg, covg_estimate, 0.9, 0.25);
}

template <class PredFn>
[[nodiscard]] static auto CallSites(std::vector<CoverageSignals> const covg,
                                    std::uint32_t const covg_estimate,
                                    PredFn pred_fn)
    -> std::enable_if_t<
        std::is_invocable_r_v<bool, PredFn, CoverageSignals const,
                              std::uint32_t const>,
        std::vector<std::uint32_t>> {
  auto buff = std::vector<std::uint32_t>();
  buff.reserve(buff.size() / 100);
  for (auto i = 0U; i < covg.size(); ++i) {
    if (pred_fn(covg[i], covg_estimate)) {
      buff.push_back(i);
    }
  }

  return decltype(buff)(buff.begin(), buff.end());
}

[[nodiscard]] static auto CallSnpCandidates(
    std::vector<CoverageSignals> const covg, std::uint32_t const covg_estimate)
    -> std::vector<std::uint32_t> {
  return CallSites(covg, covg_estimate, IsSnpSite);
}

[[nodiscard]] static auto CallErrorSites(
    std::vector<CoverageSignals> const covg,
    std::uint32_t const covg_estimate) {
  return CallSites(covg, covg_estimate, IsErrorSite);
}

[[nodiscard]] static auto CalcWindowIntervals(
    std::vector<std::uint32_t> error_sites) -> std::vector<Interval> {
  auto dst = std::vector<Interval>(error_sites.size());
  std::transform(error_sites.cbegin(), error_sites.cend(), dst.begin(),
                 [](std::uint32_t const pos) -> Interval {
                   return {pos, pos};
                 });

  auto constexpr kSentinelDist = std::numeric_limits<std::uint32_t>::max();
  auto constexpr kDeadInterval = Interval{1, 0};

  for (auto can_terminate = false; !can_terminate;) {
    can_terminate = true;
    for (auto i = 0; i < dst.size(); ++i) {
      auto const lhs_dist =
          i > 0 ? dst[i].last - dst[i - 1].first : kSentinelDist;
      auto const rhs_dist =
          i + 1 < dst.size() ? dst[i + 1].last - dst[i].first : kSentinelDist;

      if (lhs_dist != kSentinelDist || rhs_dist != kSentinelDist) {
        if (lhs_dist < rhs_dist) {
          dst[i].first = dst[i - 1].first;
          dst[i - 1] = kDeadInterval;
          can_terminate = false;
        } else {
          dst[i].last = dst[i + 1].last;
          std::swap(dst[i], dst[i + 1]);
          dst[i] = kDeadInterval;
          can_terminate = false;
        }
      }
    }

    dst.erase(std::stable_partition(dst.begin(), dst.end(),
                                    [kDeadInterval](Interval interval) -> bool {
                                      return interval != kDeadInterval;
                                    }),
              dst.end());
  }

  return dst;
}

[[nodiscard]] static auto WindowIntervalsToWindows(
    std::uint32_t const read_len, std::vector<Interval> intervals)
    -> std::vector<ReferenceWindow> {
  auto dst = std::vector<ReferenceWindow>(intervals.size());

  auto constexpr kSentinelDist = std::numeric_limits<std::uint32_t>::max();
  for (auto i = 0U; i < intervals.size(); ++i) {
    auto const default_lhs_shift = 14U;
    auto const mid_point_lhs_shift =
        (i > 0U ? (intervals[i].first - intervals[i - 1].last)
                : intervals[i].first) /
        2U;

    auto const default_rhs_shift = 14U;
    auto const mid_point_rhs_shift =
        (i + 1U < intervals.size() ? intervals[i + 1].first - intervals[i].last
                                   : read_len - 1U - intervals[i].last) /
        2U;

    auto const lhs_shift = std::min(default_lhs_shift, mid_point_lhs_shift);
    auto const rhs_shift = std::min(default_rhs_shift, mid_point_rhs_shift);

    dst[i].interval = {.first = intervals[i].first - lhs_shift,
                       .last = intervals[i].last + rhs_shift};
  }

  return dst;
}

static auto BindReadSegmentsToWindows(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads,
    std::vector<biosoup::Overlap> const& overlaps,
    std::vector<EdlibAlignResult> const& edlib_results,
    std::vector<ReferenceWindow>& windows) -> void {
  auto const cmp = [](std::uint32_t const ovlp_start,
                      ReferenceWindow const& ref_window) -> bool {
    return ovlp_start < ref_window.interval.last;
  };

  for (auto i = 0U; i < overlaps.size(); ++i) {
    auto win_idx = std::distance(
        windows.begin(),
        std::upper_bound(windows.begin(), windows.end(), overlaps[i].lhs_begin,
                         [](std::uint32_t const ovlp_start,
                            ReferenceWindow const& ref_window) -> bool {
                           return ovlp_start < ref_window.interval.last;
                         }));

    if (win_idx == windows.size()) {
      continue;
    }

    auto lhs_first = overlaps[i].lhs_begin;
    auto rhs_first = overlaps[i].rhs_begin;

    auto lhs_last = lhs_first;
    auto rhs_last = rhs_first;

    auto const rhs_view =
        NucleicView(reads[overlaps[i].rhs_id].get(), !overlaps[i].strand);
    for (auto j = 0U; j < edlib_results[i].alignmentLength; ++j) {
      if (lhs_last == windows[win_idx].interval.first) {
        lhs_first = lhs_last;
        rhs_first = rhs_last;
      } else if (lhs_last == windows[win_idx].interval.last) {
        auto aligned_interval = Interval{.first = rhs_first, .last = rhs_last};
        if (IntervalLength(aligned_interval) >= 42U) {
          windows[win_idx].aligned_segments.push_back(AlignedSegment{
              .aligned_interval = aligned_interval,
              .bases = rhs_view.InflateData(rhs_first, rhs_last - rhs_first)});
        }

        ++win_idx;
      }

      lhs_last += (edlib_results[i].alignment[j] != 2);
      rhs_last += (edlib_results[i].alignment[j] != 1);
    }
  }
}

[[nodiscard]] static auto CreateWindowsFromAlignments(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads,
    std::vector<biosoup::Overlap> const& overlaps,
    std::vector<EdlibAlignResult> const& edlib_results,
    std::uint32_t const global_coverage_estimate)
    -> std::vector<ReferenceWindow> {
  auto const query_id = overlaps.front().lhs_id;
  auto coverage = CalculateCoverage(reads, overlaps, edlib_results);
  auto windows = WindowIntervalsToWindows(
      reads[query_id]->inflated_len,
      CalcWindowIntervals(CallErrorSites(coverage, global_coverage_estimate)));
  BindReadSegmentsToWindows(reads, overlaps, edlib_results, windows);

  return windows;
}

}  // namespace detail

auto ErrorCorrect(State& state, MapCfg const map_cfg,
                  CorrectConfig const correct_cfg,
                  std::vector<std::unique_ptr<biosoup::NucleicAcid>> src_reads)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
  auto corrected_targets = std::vector<std::unique_ptr<biosoup::NucleicAcid>>();
  auto dst = std::vector<std::unique_ptr<biosoup::NucleicAcid>>();
  auto timer = biosoup::Timer();

  timer.Start();
  auto const overlaps = FindOverlaps(state, map_cfg, src_reads);

  auto target_ids = std::vector<std::uint32_t>();
  for (auto read_id = 0U; read_id < overlaps.size(); ++read_id) {
    if (!overlaps[read_id].empty()) {
      target_ids.push_back(read_id);
    }
  }

  auto const kNTargets = target_ids.size();
  auto const kNContained = overlaps.size() - kNTargets;

  auto const kNOverlaps = std::transform_reduce(
      overlaps.cbegin(), overlaps.cend(), 0UL, std::plus<std::size_t>(),
      [](std::vector<biosoup::Overlap> const& ovlps) -> std::size_t {
        return ovlps.size();
      });

  timer.Stop();

  timer.Start();
  auto const kCovgEstimate =
      detail::EstimateCoverage(state, src_reads, overlaps);
  fmt::print(stderr, "[camel::ErrorCorrect]({:12.3f}) coverage estimate: {}\n",
             timer.Stop(), kCovgEstimate);

  {
    timer.Start();
    auto align_futures =
        std::vector<std::vector<std::future<EdlibAlignResult>>>(kNTargets);

    for (auto i = 0U; i < kNTargets; ++i) {
      auto const target_id = target_ids[i];
      align_futures[i].resize(overlaps[target_id].size());
      for (auto j = 0U; j < overlaps[target_id].size(); ++j) {
        align_futures[i][j] = state.thread_pool->Submit(
            [&src_reads, &ovlp = overlaps[target_id][j]]() -> EdlibAlignResult {
              auto const [lhs_str, rhs_str] =
                  detail::ExtractSubstrings(src_reads, ovlp);
              return detail::AlignStrings(lhs_str, rhs_str);
            });
      }
    }

    auto window_futures =
        std::vector<std::future<std::vector<detail::ReferenceWindow>>>(
            align_futures.size());

    auto n_transformed = 0U;
    for (auto i = 0U; i < kNTargets; ++i) {
      auto alignments = std::vector<EdlibAlignResult>(align_futures[i].size());
      std::transform(std::make_move_iterator(align_futures[i].begin()),
                     std::make_move_iterator(align_futures[i].end()),
                     alignments.begin(),
                     std::mem_fn(&std::future<EdlibAlignResult>::get));

      n_transformed += alignments.size();
      window_futures[i] = state.thread_pool->Submit(
          [&src_reads, &overlaps, kCovgEstimate,
           alignments = std::move(alignments)](
              std::uint32_t read_id) -> std::vector<detail::ReferenceWindow> {
            auto windows = detail::CreateWindowsFromAlignments(
                src_reads, overlaps[read_id], alignments, kCovgEstimate);

            std::for_each(alignments.begin(), alignments.end(),
                          edlibFreeAlignResult);

            return windows;
          },
          target_ids[i]);

      fmt::print(stderr,
                 "\r[camel::ErrorCorrect({:12.3f})] transformed {} / {} "
                 "alignment futures to window tasks",
                 timer.Lap(), n_transformed, kNOverlaps);
    }

    decltype(align_futures){}.swap(align_futures);
    fmt::print(stderr,
               "\r[camel::ErrorCorrect]({:12.3f}) transformed {} / {} "
               "alignment futures "
               "to window tasks\n",
               timer.Lap(), n_transformed, kNOverlaps);

    for (auto i = 0U; i < window_futures.size(); ++i) {
      window_futures[i].wait();
      if ((i & 127) == 0U) {
        fmt::print(stderr,
                   "\r[camel::ErrorCorrect]({:12.3f}) collected window futures "
                   "{} / {}",
                   timer.Lap(), i + 1, window_futures.size());
      }
    }

    fmt::print(stderr,
               "\r[camel::ErrorCorrect]({:12.3f}) collected window futures "
               "{} / {}\n",
               timer.Lap(), window_futures.size(), window_futures.size());
  }

  fmt::print(stderr, "[camel::ErrorCorrect]({:12.3f})\n", timer.elapsed_time());
  return dst;
}

}  // namespace camel
