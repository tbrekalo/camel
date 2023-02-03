#include "window.h"

#include <algorithm>
#include <limits>

#include "call_sites.h"
#include "spoa/alignment_engine.hpp"
#include "spoa/graph.hpp"
#include "tbb/enumerable_thread_specific.h"

namespace camel::detail {

static constexpr auto kInvalidIndex = std::numeric_limits<uint32_t>::max();

static auto AlignmentEngines =
    tbb::enumerable_thread_specific<std::unique_ptr<spoa::AlignmentEngine>>();

static auto GetAlignmentEngine(POAConfig const config)
    -> spoa::AlignmentEngine& {
  auto& engine = AlignmentEngines.local();
  if (!engine) {
    engine = spoa::AlignmentEngine::Create(
        spoa::AlignmentType::kNW, config.match, config.mismatch, config.gap);
  }

  return *engine;
}

static auto BindReadSegmentsToWindows(
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> reads,
    std::span<biosoup::Overlap const> overlaps,
    std::span<EdlibAlignResult const> edlib_results,
    std::span<ReferenceWindow> windows) -> void {
  for (auto ovlp_idx = 0U; ovlp_idx < overlaps.size(); ++ovlp_idx) {
    auto win_idx = std::distance(
        windows.begin(),
        std::upper_bound(windows.begin(), windows.end(),
                         overlaps[ovlp_idx].rhs_begin,
                         [](std::uint32_t const ovlp_start,
                            ReferenceWindow const& ref_window) -> bool {
                           return ovlp_start < ref_window.interval.last;
                         }));

    if (win_idx >= windows.size()) {
      continue;
    }

    auto const query_view = NucleicView(reads[overlaps[ovlp_idx].lhs_id].get(),
                                        !overlaps[ovlp_idx].strand);

    auto query_first = kInvalidIndex;
    auto target_first = kInvalidIndex;

    auto query_last = query_first;
    auto target_last = target_first;
    auto query_curr = overlaps[ovlp_idx].strand
                          ? overlaps[ovlp_idx].lhs_begin
                          : reads[overlaps[ovlp_idx].lhs_id]->inflated_len -
                                overlaps[ovlp_idx].lhs_end;
    auto target_curr = overlaps[ovlp_idx].rhs_begin;

    for (auto align_idx = 0U;
         align_idx < edlib_results[ovlp_idx].alignmentLength; ++align_idx) {
      if (target_curr >= windows[win_idx].interval.first &&
          (edlib_results[ovlp_idx].alignment[align_idx] == 0 ||
           edlib_results[ovlp_idx].alignment[align_idx] == 3)) {
        if (query_first == kInvalidIndex) {
          query_first = query_curr;
          target_first = target_curr;
        }

        query_last = query_curr + 1;
        target_last = target_curr + 1;
      }

      query_curr += (edlib_results[ovlp_idx].alignment[align_idx] != 2);
      target_curr += (edlib_results[ovlp_idx].alignment[align_idx] != 1);

      if (target_curr != windows[win_idx].interval.last &&
          align_idx + 1 != edlib_results[ovlp_idx].alignmentLength) {
        continue;
      }

      auto quality_sum = 0.0;
      for (auto pos = query_first; pos < query_curr; ++pos) {
        quality_sum += query_view.Quality(pos) - 33;
      }

      if (quality_sum / (query_last - query_first) > 10) {
        windows[win_idx].aligned_segments.emplace_back(AlignedSegment{
            .alignment_local_interval = LocalizeInterval(
                windows[win_idx].interval.first, {target_first, target_last}),
            .bases =
                query_view.InflateData(query_first, query_last - query_first),
            .quality = query_view.InflateQuality(query_first,
                                                 query_last - query_first)});
      }

      query_first = target_first = kInvalidIndex;
      if (++win_idx >= windows.size()) {
        break;
      }
    }
  }

  for (auto& win : windows) {
    std::sort(win.aligned_segments.begin(), win.aligned_segments.end(),
              [](AlignedSegment const& lhs, AlignedSegment const& rhs) -> bool {
                return lhs.alignment_local_interval.first <
                       rhs.alignment_local_interval.first;
              });
  }
}

static auto KeepHaploidOverlaps(
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> reads,
    std::span<biosoup::Overlap const> overlaps,
    std::span<EdlibAlignResult const> edlib_results,
    std::span<CoverageSignals const> coverage,
    std::uint32_t const global_coverage_estimate)
    -> std::pair<std::vector<biosoup::Overlap>, std::vector<EdlibAlignResult>> {
  auto haploid_overlaps = std::vector<biosoup::Overlap>();
  auto haploid_alignments = std::vector<EdlibAlignResult>();

  auto const target_id = overlaps.front().rhs_id;
  auto snp_sites = std::vector<std::uint32_t>();
  snp_sites.reserve(coverage.size() / 10U);

  auto read_view = NucleicView(reads[target_id].get(), false);

  for (auto i = 0U; i < read_view.InflatedLenght(); ++i) {
    if (IsStableSite(coverage[i], global_coverage_estimate, read_view.Code(i),
                     0.25, 0.5) &&
        IsSnpSite(coverage[i], global_coverage_estimate, 0.8)) {
      snp_sites.push_back(i);
    }
  }

  if (!snp_sites.empty()) {
    auto normed_snp_rate = std::vector<double>(overlaps.size(), 0);
    for (auto ovlp_idx = 0U; ovlp_idx < overlaps.size(); ++ovlp_idx) {
      auto const& edlib_res = edlib_results[ovlp_idx];
      auto target_pos = overlaps[ovlp_idx].rhs_begin;

      auto snp_site_idx = 0U;
      while (snp_site_idx < snp_sites.size() &&
             snp_sites[snp_site_idx] < target_pos) {
        ++snp_site_idx;
      }
      for (auto i = 0;
           i < edlib_res.alignmentLength && snp_site_idx < snp_sites.size();
           ++i) {
        target_pos += (edlib_res.alignment[i] != 1);
        while (snp_site_idx < snp_sites.size() &&
               snp_sites[snp_site_idx] < target_pos) {
          ++snp_site_idx;
        }
        if (snp_site_idx < snp_sites.size() &&
            snp_sites[snp_site_idx] == target_pos) {
          normed_snp_rate[ovlp_idx] += edlib_res.alignment[i] != 0;
        }
      }
    }

    for (auto i = 0U; i < overlaps.size(); ++i) {
      normed_snp_rate[i] /= overlaps[i].rhs_end - overlaps[i].rhs_begin;
      normed_snp_rate[i] *= reads[overlaps[i].rhs_id]->inflated_len;
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
          return lhs.first < rhs.first;
        });

    auto const median_snp_rate =
        scores_indices[scores_indices.size() / 2].first;

    auto const n_kept_ovlps =
        std::count_if(normed_snp_rate.cbegin(), normed_snp_rate.cend(),
                      [median_snp_rate](double const snp_rate) -> bool {
                        return snp_rate <= median_snp_rate;
                      });

    haploid_overlaps.resize(n_kept_ovlps);
    haploid_alignments.resize(n_kept_ovlps);

    auto j = 0U;
    for (auto i = 0U; i < overlaps.size(); ++i) {
      if (j < n_kept_ovlps && normed_snp_rate[i] <= median_snp_rate) {
        haploid_overlaps[j] = overlaps[i];
        haploid_alignments[j++] = edlib_results[i];
      } else {
        edlibFreeAlignResult(edlib_results[i]);
      }
    }
  } else {
    haploid_overlaps.insert(haploid_overlaps.end(), overlaps.begin(),
                            overlaps.end());
    haploid_alignments.insert(haploid_alignments.end(), edlib_results.begin(),
                              edlib_results.end());
  }

  return std::make_pair(std::move(haploid_overlaps),
                        std::move(haploid_alignments));
}

auto FindWindows(NucleicView read_view, std::span<CoverageSignals> coverage,
                 std::uint32_t const global_coverage_estimate,
                 std::uint32_t const window_len)
    -> std::vector<ReferenceWindow> {
  auto dst = std::vector<ReferenceWindow>();

  for (auto i = 0U; i < read_view.InflatedLenght();) {
    auto j = std::min(i + window_len, read_view.InflatedLenght());
    dst.push_back(ReferenceWindow{.interval = {i, j}});
    i = j;
  }

  return dst;
}

auto CreateWindowsFromAlignments(
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> reads,
    std::span<biosoup::Overlap const> overlaps,
    std::span<EdlibAlignResult const> alignments,
    std::uint32_t const window_len,
    std::uint32_t const global_coverage_estimate)
    -> std::vector<ReferenceWindow> {
  auto const target_id = overlaps.front().rhs_id;
  auto coverage = CalculateCoverage(reads, overlaps, alignments);

  auto [haploid_overlaps, haploid_alignments] = KeepHaploidOverlaps(
      reads, overlaps, alignments, coverage, global_coverage_estimate);

  auto read_view = NucleicView(reads[target_id].get(), false);
  auto windows =
      FindWindows(read_view, coverage, global_coverage_estimate, window_len);

  BindReadSegmentsToWindows(reads, haploid_overlaps, haploid_alignments,
                            windows);

  std::for_each(haploid_alignments.begin(), haploid_alignments.end(),
                edlibFreeAlignResult);
  return windows;
}

auto ReleaseAlignmentEngines() -> std::size_t {
  auto const dst = AlignmentEngines.size();
  AlignmentEngines.clear();

  return dst;
}

auto WindowConsensus(std::string_view backbone_data,
                     std::string_view backbone_quality,
                     ReferenceWindowView ref_window_view, POAConfig poa_cfg)
    -> ConsensusResult {
  if (ref_window_view.aligned_segments.size() < 3) {
    return {.bases = std::string(backbone_data.begin(), backbone_data.end()),
            .is_corrected = false};
  }

  auto& alignment_engine = detail::GetAlignmentEngine(poa_cfg);
  auto graph = spoa::Graph();

  {
    auto local_backbone_str = std::string(backbone_data);
    auto local_backbone_quality_str = std::string(backbone_quality);
    if (local_backbone_quality_str.empty()) {
      graph.AddAlignment(spoa::Alignment(), local_backbone_str);
    } else {
      graph.AddAlignment(spoa::Alignment(), local_backbone_str,
                         local_backbone_quality_str);
    }
  }

  auto const [win_ref_intv, aligned_segments] = ref_window_view;
  auto const window_length = IntervalLength(win_ref_intv);

  auto const legal_start = static_cast<decltype(window_length)>(
      window_length * detail::kAllowedFuzzPercent);
  auto const legal_end = window_length - legal_start;
  for (auto const& [alignment_interval, bases, quality] : aligned_segments) {
    auto const alignment_len = IntervalLength(alignment_interval);
    if (alignment_len < window_length * detail::kSmallWindowPercent) {
      continue;
    }

    auto alignment = spoa::Alignment();
    if (alignment_interval.first < legal_start &&
        alignment_interval.last - 1 > legal_end) {
      alignment = alignment_engine.Align(bases, graph);

    } else {
      auto mapping = std::vector<spoa::Graph::Node const*>();
      auto subgraph = graph.Subgraph(alignment_interval.first,
                                     alignment_interval.last - 1, &mapping);

      alignment = alignment_engine.Align(bases, subgraph);
      subgraph.UpdateAlignment(mapping, &alignment);
    }

    graph.AddAlignment(alignment, bases, quality);
  }

  std::vector<std::uint32_t> coverages;
  auto consensus = graph.GenerateConsensus(&coverages);
  uint32_t average_coverage = (ref_window_view.aligned_segments.size() - 1) / 2;

  int32_t begin = 0, end = consensus.size() - 1;
  for (; begin < static_cast<int32_t>(consensus.size()); ++begin) {
    if (coverages[begin] >= average_coverage) {
      break;
    }
  }
  for (; end >= 0; --end) {
    if (coverages[end] >= average_coverage) {
      break;
    }
  }

  if (begin < end) {
    consensus = consensus.substr(begin, end - begin + 1);
  }

  return {.bases = std::move(consensus), .is_corrected = true};
}

}  // namespace camel::detail
