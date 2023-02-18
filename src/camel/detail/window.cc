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
    std::span<ReferenceWindow> windows, double quality_threshold) -> void {
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

    auto const& cigar = overlaps[ovlp_idx].alignment;
    auto const query_view = NucleicView(reads[overlaps[ovlp_idx].lhs_id].get(),
                                        !overlaps[ovlp_idx].strand);

    auto query_first = kInvalidIndex, target_first = kInvalidIndex;
    auto query_last = query_first, target_last = target_first;

    auto query_curr = (overlaps[ovlp_idx].strand
                           ? overlaps[ovlp_idx].lhs_begin
                           : reads[overlaps[ovlp_idx].lhs_id]->inflated_len -
                                 overlaps[ovlp_idx].lhs_end) -
                      1;
    auto target_curr = overlaps[ovlp_idx].rhs_begin - 1;

    auto const update_window = [&]() {
      if (target_first != kInvalidIndex) {
        auto quality_sum = 0.0;
        for (auto pos = query_first; pos < query_last; ++pos) {
          quality_sum += query_view.Quality(pos) - 33;
        }

        if (quality_sum / (query_last - query_first) > quality_threshold) {
          windows[win_idx].aligned_segments.emplace_back(AlignedSegment{
              .alignment_local_interval = LocalizeInterval(
                  windows[win_idx].interval.first, {target_first, target_last}),
              .bases =
                  query_view.InflateData(query_first, query_last - query_first),
              .quality = query_view.InflateQuality(query_first,
                                                   query_last - query_first)});
        }
      }

      query_first = target_first = kInvalidIndex;
      return ++win_idx == windows.size();
    };

    auto windowsTraversed = false;
    for (auto j = 0U, i = 0U; i < cigar.size() && !windowsTraversed; ++i) {
      if (IsMisOrMatch(cigar[i])) {
        auto const n = std::atoi(&cigar[j]);
        for (auto k = 0; k < n && !windowsTraversed; ++k) {
          ++query_curr;
          ++target_curr;
          if (target_first == kInvalidIndex) {
            query_first = query_curr;
            target_first = target_curr;
          }

          query_last = query_curr + 1;
          target_last = target_curr + 1;

          if (target_curr + 1 == windows[win_idx].interval.last) {
            windowsTraversed = update_window();
          }
        }
        j = i + 1;
      } else if (IsInsertion(cigar[i])) {
        query_curr += std::atoi(&cigar[j]);
        j = i + 1;
      } else if (IsDeletion(cigar[i])) {
        auto const n = std::atoi(&cigar[j]);
        for (auto k = 0; k < n && !windowsTraversed; ++k) {
          ++target_curr;
          if (target_first != kInvalidIndex) {
            target_last = target_curr;
            if (target_curr == windows[win_idx].interval.last) {
              windowsTraversed = update_window();
            }
          }
        }
        j = i + 1;
      } else if (IsClipOrPad(cigar[i])) {
        j = i + 1;
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

auto FindWindows(NucleicView read_view,
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
    std::span<biosoup::Overlap const> overlaps, WindowConfig window_cfg,
    std::uint32_t const global_coverage_estimate)
    -> std::vector<ReferenceWindow> {
  auto const target_id = overlaps.front().rhs_id;

  auto read_view = NucleicView(reads[target_id].get(), false);
  auto windows = FindWindows(read_view, global_coverage_estimate,
                             window_cfg.window_length);

  BindReadSegmentsToWindows(reads, overlaps, windows,
                            window_cfg.quality_threshold);

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
  uint32_t average_coverage = ref_window_view.aligned_segments.size() / 2;

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
