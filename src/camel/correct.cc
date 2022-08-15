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
static auto constexpr kWinLength = 240U;

static auto constexpr kWinPadding = 13U;
static auto constexpr kAllowedFuzzPercent = 0.05;
static auto constexpr kSmallWindowPercent = 0.05;

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

constexpr auto LocalizeInterval(std::uint32_t const pos, Interval const intv)
    -> Interval {
  return {.first = intv.first - pos, .last = intv.last - pos};
}

struct AlignedSegment {
  Interval alignment_local_interval;
  std::string bases;
};

struct ReferenceWindow {
  Interval interval;
  std::vector<AlignedSegment> aligned_segments;
};

auto operator<<(std::ostream& ostrm, ReferenceWindow ref_win) -> std::ostream& {
  return ostrm << IntervalLength(ref_win.interval);
}

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

    auto i = 0;
    for (; i < edlib_res.alignmentLength && edlib_res.alignment[i] == 2; ++i)
      ;  // skip initial insertion
    for (; i < edlib_res.alignmentLength; ++i) {
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
          ++dst[lhs_pos - 1U].signals[CoverageSignals::kInsIdx];
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

  if (*std::max_element(covg.signals.cbegin(), covg.signals.cend()) <
      covg_estimate / 3) {
    return false;
  }

  auto const strong_base_threshold =
      static_cast<std::uint32_t>(std::round(strong_base_ratio * sig_sum));
  auto const indle_threshold =
      static_cast<std::uint32_t>(std::round(indle_ratio * sig_sum));
  auto const cutoff_threshold =
      static_cast<std::uint32_t>(std::round(covg_estimate * kSqrt5));

  auto dst = true;
  for (auto i = 0U; dst && i < 4U; ++i) {
    if (covg.signals[i] >= strong_base_threshold) {
      dst = false;
    }
  }

  for (auto i = 4U; dst && i < 6U; ++i) {
    if (covg.signals[i] > indle_threshold) {
      dst = false;
    }
  }

  return dst;
}

[[nodiscard]] static auto IsSnpSite(CoverageSignals const& covg,
                                    std::uint32_t const covg_estimate) -> bool {
  auto signals = covg.signals;
  std::sort(signals.begin(), std::next(signals.begin(), 4),
            std::greater<std::uint16_t>());

  auto const sig_sum = std::accumulate(signals.cbegin(), signals.cend(), 0U,
                                       std::plus<std::uint32_t>());

  if (sig_sum > 0.5 * covg_estimate && signals[0] < 0.8 * sig_sum &&
      signals[1] > 0.2 * sig_sum && signals[0] + signals[1] >= 0.9 * sig_sum) {
    return true;
  }

  return false;
}

[[nodiscard]] static auto IsErrorSite(CoverageSignals const& covg,
                                      std::uint32_t const covg_estimate)
    -> bool {
  return IsUnstableSite(covg, covg_estimate, 0.95, 0.25);
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
      auto const lhs_dist = i > 0 && dst[i - 1] != kDeadInterval
                                ? dst[i].last - dst[i - 1].first
                                : kSentinelDist;
      auto const rhs_dist = i + 1 < dst.size() && dst[i + 1] != kDeadInterval
                                ? dst[i + 1].last - dst[i].first
                                : kSentinelDist;

      if (lhs_dist <= kWinLength || rhs_dist <= kWinLength) {
        if (lhs_dist < rhs_dist) {
          dst[i].first = dst[i - 1].first;
          dst[i - 1] = kDeadInterval;
        } else {
          dst[i].last = dst[i + 1].last;
          std::swap(dst[i], dst[i + 1]);

          dst[i] = kDeadInterval;
          if (i > 0) {
            std::swap(dst[i - 1], dst[i]);
          }
        }
        can_terminate = false;
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
  for (auto i = 0U; i < overlaps.size(); ++i) {
    auto win_idx = std::distance(
        windows.begin(),
        std::upper_bound(windows.begin(), windows.end(), overlaps[i].lhs_begin,
                         [](std::uint32_t const ovlp_start,
                            ReferenceWindow const& ref_window) -> bool {
                           return ovlp_start < ref_window.interval.last;
                         }));

    if (win_idx >= windows.size()) {
      continue;
    }

    auto lhs_first = overlaps[i].lhs_begin;
    auto rhs_first = overlaps[i].rhs_begin;

    auto lhs_curr = lhs_first;
    auto rhs_curr = rhs_first;

    auto const rhs_view =
        NucleicView(reads[overlaps[i].rhs_id].get(), !overlaps[i].strand);
    for (auto j = 0U; j < edlib_results[i].alignmentLength; ++j) {
      if (lhs_curr == windows[win_idx].interval.first) {
        lhs_first = lhs_curr;
        rhs_first = rhs_curr;
      }

      lhs_curr += (edlib_results[i].alignment[j] != 2);
      rhs_curr += (edlib_results[i].alignment[j] != 1);

      if (lhs_curr == windows[win_idx].interval.last) {
        windows[win_idx].aligned_segments.emplace_back(AlignedSegment{
            .alignment_local_interval = LocalizeInterval(
                windows[win_idx].interval.first, {lhs_first, lhs_curr}),
            .bases = rhs_view.InflateData(rhs_first, rhs_curr - rhs_first)});

        if (++win_idx >= windows.size()) {
          break;
        }
      }
    }
  }
}

static auto KeepHaploidOverlaps(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads,
    std::vector<biosoup::Overlap>& overlaps,
    std::vector<EdlibAlignResult>& edlib_results,
    std::vector<CoverageSignals> const& coverage,
    std::uint32_t const global_coverage_estimate) -> void {
  auto const query_id = overlaps.front().lhs_id;
  auto snp_sites = std::vector<std::uint32_t>();
  snp_sites.reserve(coverage.size() / 10U);

  for (auto i = 0U; i < coverage.size(); ++i) {
    if (IsSnpSite(coverage[i], global_coverage_estimate)) {
      snp_sites.push_back(i);
    }
  }

  if (!snp_sites.empty()) {
    auto normed_snp_rate = std::vector<double>(overlaps.size(), 0);
    for (auto ovlp_idx = 0U; ovlp_idx < overlaps.size(); ++ovlp_idx) {
      auto const& edlib_res = edlib_results[ovlp_idx];

      auto lhs_pos = overlaps[ovlp_idx].lhs_begin;

      auto i = 0;
      auto snp_site_idx = 0U;
      for (; i < edlib_res.alignmentLength && edlib_res.alignment[i] == 2; ++i)
        ;  // skip initial insertion
      while (snp_site_idx < snp_sites.size() &&
             snp_sites[snp_site_idx] < lhs_pos) {
        ++snp_site_idx;
      }
      for (; i < edlib_res.alignmentLength && snp_site_idx < snp_sites.size();
           ++i) {
        lhs_pos += (edlib_res.alignment[i] != 2);
        while (snp_site_idx < snp_sites.size() &&
               snp_sites[snp_site_idx] < lhs_pos) {
          ++snp_site_idx;
        }
        if (snp_sites[snp_site_idx] == lhs_pos) {
          normed_snp_rate[ovlp_idx] += edlib_res.alignment[i] != 0;
        }
      }
    }

    for (auto i = 0U; i < overlaps.size(); ++i) {
      normed_snp_rate[i] /= overlaps[i].rhs_end - overlaps[i].rhs_begin;
      normed_snp_rate[i] *= reads[overlaps[i].lhs_id]->inflated_len;
    }

    auto scores_indices =
        std::vector<std::pair<double, std::uint32_t>>(normed_snp_rate.size());

    for (auto i = 0U; i < normed_snp_rate.size(); ++i) {
      scores_indices[i] = {normed_snp_rate[i], i};
    }

    std::nth_element(
        scores_indices.begin(),
        std::next(scores_indices.begin(), scores_indices.size() / 2),
        scores_indices.end(),
        [](std::pair<double, std::uint32_t> const lhs,
           std::pair<double, std::uint32_t> const rhs) -> bool {
          return lhs.second < rhs.second;
        });

    auto const median_snp_rate =
        scores_indices[scores_indices.size() / 2].first;

    auto const n_kept_ovlps =
        std::count_if(normed_snp_rate.cbegin(), normed_snp_rate.cend(),
                      [median_snp_rate](double const snp_rate) -> bool {
                        return snp_rate <= median_snp_rate;
                      });

    auto updated_ovlps = std::vector<biosoup::Overlap>(n_kept_ovlps);
    auto updated_edlib_results = std::vector<EdlibAlignResult>(n_kept_ovlps);

    auto j = 0U;
    for (auto i = 0U; i < overlaps.size(); ++i) {
      if (j < n_kept_ovlps && normed_snp_rate[i] <= median_snp_rate) {
        updated_ovlps[j] = overlaps[i];
        updated_edlib_results[j++] = edlib_results[i];
      } else {
        edlibFreeAlignResult(edlib_results[i]);
      }
    }

    std::swap(updated_ovlps, overlaps);
    std::swap(updated_edlib_results, edlib_results);
  }
}

[[nodiscard]] static auto CreateWindowsFromAlignments(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads,
    std::vector<biosoup::Overlap> overlaps,
    std::vector<EdlibAlignResult> edlib_results,
    std::uint32_t const global_coverage_estimate)
    -> std::vector<ReferenceWindow> {
  auto const query_id = overlaps.front().lhs_id;
  auto coverage = CalculateCoverage(reads, overlaps, edlib_results);

  KeepHaploidOverlaps(reads, overlaps, edlib_results, coverage,
                      global_coverage_estimate);

  auto windows = WindowIntervalsToWindows(
      reads[query_id]->inflated_len,
      CalcWindowIntervals(CallErrorSites(coverage, global_coverage_estimate)));

  BindReadSegmentsToWindows(reads, overlaps, edlib_results, windows);

  std::for_each(edlib_results.begin(), edlib_results.end(),
                edlibFreeAlignResult);

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
                src_reads, std::move(overlaps[read_id]), std::move(alignments),
                kCovgEstimate);

            return windows;
          },
          target_ids[i]);

      fmt::print(stderr,
                 "\r[camel::ErrorCorrect({:12.3f})] transformed {} / {} "
                 "alignment futures to window tasks",
                 timer.Lap(), n_transformed, kNOverlaps);
    }

    fmt::print(stderr,
               "\r[camel::ErrorCorrect]({:12.3f}) transformed {} / {} "
               "alignment futures "
               "to window tasks\n",
               timer.Lap(), n_transformed, kNOverlaps);
    decltype(align_futures){}.swap(align_futures);

    auto ref_windows = std::vector<std::vector<detail::ReferenceWindow>>();
    ref_windows.reserve(window_futures.size());

    for (auto i = 0U; i < window_futures.size(); ++i) {
      ref_windows.push_back(window_futures[i].get());

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
               timer.Stop(), window_futures.size(), window_futures.size());
    decltype(window_futures)().swap(window_futures);

    auto alignment_engines =
        tsl::robin_map<std::thread::id,
                       std::unique_ptr<spoa::AlignmentEngine>>();

    for (auto const [thread_id, _] : state.thread_pool->thread_map()) {
      alignment_engines[thread_id] = spoa::AlignmentEngine::Create(
          spoa::AlignmentType::kNW, correct_cfg.poa_cfg.match,
          correct_cfg.poa_cfg.mismatch, correct_cfg.poa_cfg.gap);
    }

    timer.Start();
    dst.reserve(window_futures.size());
    auto spoa_futures = std::vector<std::future<void>>();

    for (auto read_idx = 0U; read_idx < ref_windows.size(); ++read_idx) {
      auto graphs = std::vector<spoa::Graph>(ref_windows[read_idx].size());
      auto backbone = src_reads[target_ids[read_idx]]->InflateData();

      for (auto win_idx = 0U; win_idx < ref_windows[read_idx].size();
           ++win_idx) {
        graphs[win_idx].AddAlignment(
            spoa::Alignment(),
            backbone.substr(ref_windows[read_idx][win_idx].interval.first,
                            detail::IntervalLength(
                                ref_windows[read_idx][win_idx].interval)));
      }

      for (auto win_idx = 0U; win_idx < ref_windows[read_idx].size();
           ++win_idx) {
        spoa_futures.emplace_back(state.thread_pool->Submit(
            [&alignment_engines,
             active_window = std::move(ref_windows[read_idx][win_idx]), &graphs,
             win_idx]() -> void {
              auto const& alignment_engine =
                  alignment_engines[std::this_thread::get_id()];

              auto const& [win_ref_intv, aligned_segments] = active_window;
              auto const window_length = IntervalLength(win_ref_intv);

              for (auto const& [alignment_interval, bases] : aligned_segments) {
                auto const alignment_len = IntervalLength(alignment_interval);
                if (alignment_len <
                    window_length * detail::kSmallWindowPercent) {
                  continue;
                }

                auto const legal_start =
                    window_length * detail::kAllowedFuzzPercent;
                auto const legal_end = window_length - legal_start;

                auto alignment = spoa::Alignment();
                if (alignment_interval.first <= legal_start &&
                    legal_end <= alignment_interval.last) {
                  alignment = alignment_engine->Align(bases, graphs[win_idx]);
                } else {
                  auto mapping = std::vector<spoa::Graph::Node const*>();
                  auto subgraph = graphs[win_idx].Subgraph(
                      alignment_interval.first, alignment_interval.last - 1,
                      &mapping);

                  alignment =
                      alignment_engines[std::this_thread::get_id()]->Align(
                          bases, subgraph);
                  subgraph.UpdateAlignment(mapping, &alignment);
                }

                graphs[win_idx].AddAlignment(alignment, bases);
              }
            }));
      }

      auto consensus = std::string();
      consensus.reserve(src_reads[target_ids[read_idx]]->inflated_len * 1.2);

      for (auto i = 0; i < spoa_futures.size(); ++i) {
        spoa_futures[i].wait();
      }

      {
        auto prev = 0U;
        for (auto win_idx = 0U; win_idx < spoa_futures.size(); ++win_idx) {
          spoa_futures[win_idx].wait();
          auto const& interval = ref_windows[read_idx][win_idx].interval;

          consensus.insert(consensus.end(), std::next(backbone.begin(), prev),
                           std::next(backbone.begin(), interval.first));
          consensus += graphs[win_idx].GenerateConsensus();
          prev = interval.last;
        }

        consensus.insert(consensus.end(), std::next(backbone.begin(), prev),
                         backbone.end());
      }

      dst.push_back(std::make_unique<biosoup::NucleicAcid>(
          src_reads[target_ids[read_idx]]->name, consensus));
      fmt::print(
          stderr,
          "\r[camel::ErrorCorrect]({:12.3f}) generated consensus for {} / "
          "{} reads",
          timer.Lap(), read_idx + 1, ref_windows.size());
      spoa_futures.clear();
    }
    fmt::print(stderr,
               "\r[camel::ErrorCorrect]({:12.3f}) generated consensus for {} / "
               "{} reads\n",
               timer.Stop(), ref_windows.size(), ref_windows.size());
  }

  fmt::print(stderr, "[camel::ErrorCorrect]({:12.3f})\n", timer.elapsed_time());
  return dst;
}

}  // namespace camel
