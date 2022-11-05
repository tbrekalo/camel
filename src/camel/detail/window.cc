#include "window.h"

#include <algorithm>

#include "call_sites.h"
#include "spoa/alignment_engine.hpp"
#include "spoa/graph.hpp"
#include "tbb/enumerable_thread_specific.h"

namespace camel::detail {

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
            .bases = rhs_view.InflateData(rhs_first, rhs_curr - rhs_first),
            .quality =
                rhs_view.InflateQuality(rhs_first, rhs_curr - rhs_first)});

        if (++win_idx >= windows.size()) {
          break;
        }
      }
    }
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

  auto const query_id = overlaps.front().lhs_id;
  auto snp_sites = std::vector<std::uint32_t>();
  snp_sites.reserve(coverage.size() / 10U);

  auto read_view = NucleicView(reads[query_id].get(), false);

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
        if (snp_site_idx < snp_sites.size() &&
            snp_sites[snp_site_idx] == lhs_pos) {
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
  auto const look_behind_anchor =
      [=](std::uint32_t const last) -> std::uint32_t {
    auto interval = Interval{last, last};
    for (auto first = last - 1; first > last && last - first < window_len / 2;
         --first) {
      if (IsStableSite(coverage[first], global_coverage_estimate,
                       read_view.Code(first), 0.7, 0.3)) {
        interval.first = first;
      } else {
        interval = {first, first};
      }

      if (IntervalLength(interval) > 14U) {
        return interval.first + IntervalLength(interval) / 2;
      }
    }
    return last;
  };

  for (auto i = 0U; i < read_view.InflatedLenght();) {
    auto j = std::min(look_behind_anchor(i + window_len),
                      read_view.InflatedLenght());
    dst.push_back(ReferenceWindow{.interval = {i, j}});
    i = j;

    if (read_view.InflatedLenght() - i < window_len / 5) {
      dst.back().interval.last = read_view.InflatedLenght();
      break;
    }
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
  auto const query_id = overlaps.front().lhs_id;
  auto coverage = CalculateCoverage(reads, overlaps, alignments);

  auto [haploid_overlaps, haploid_alignments] = KeepHaploidOverlaps(
      reads, overlaps, alignments, coverage, global_coverage_estimate);

  auto read_view = NucleicView(reads[query_id].get(), false);
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
    -> std::string {
  if (ref_window_view.aligned_segments.size() < 3) {
    return std::string(backbone_data.begin(), backbone_data.end());
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

  for (auto const& [alignment_interval, bases, quality] : aligned_segments) {
    auto const alignment_len = IntervalLength(alignment_interval);
    if (alignment_len < window_length * detail::kSmallWindowPercent) {
      continue;
    }

    auto const legal_start = window_length * detail::kAllowedFuzzPercent;
    auto const legal_end = window_length - legal_start;

    auto alignment = spoa::Alignment();
    if (alignment_interval.first <= legal_start &&
        legal_end <= alignment_interval.last) {
      alignment = alignment_engine.Align(bases, graph);
    } else if (IntervalLength(alignment_interval) > kWinPadding) {
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

  return consensus;
}

}  // namespace camel::detail
