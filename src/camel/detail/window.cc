#include "window.h"

#include <algorithm>
#include <limits>

#include "call_sites.h"
#include "spoa/alignment_engine.hpp"
#include "spoa/graph.hpp"
#include "tbb/enumerable_thread_specific.h"

namespace camel::detail {

static constexpr auto kInvalidIndex = std::numeric_limits<uint32_t>::max() - 1U;

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

struct ReadWindowIndices {
  std::int32_t curr;
  Interval active_interval;
};

struct WindowIndices {
  ReadWindowIndices target_indices;
  ReadWindowIndices query_indices;
};

static auto MakeWindowIndices(biosoup::Overlap const& ovlp,
                              std::uint32_t const query_len) -> WindowIndices {
  return WindowIndices{
      .target_indices = {.curr = static_cast<std::int32_t>(ovlp.rhs_begin) - 1,
                         .active_interval = {kInvalidIndex, kInvalidIndex}},
      .query_indices = {
          .curr = static_cast<std::int32_t>(
                      ovlp.strand ? ovlp.lhs_begin : query_len - ovlp.lhs_end) -
                  1,
          .active_interval = {kInvalidIndex, kInvalidIndex}}};
}

static auto IncrementQueryCurr(WindowIndices wi) -> WindowIndices {
  ++wi.query_indices.curr;
  return wi;
}
static auto IncrementQueryCurr(WindowIndices wi, std::uint32_t const val)
    -> WindowIndices {
  wi.query_indices.curr += val;
  return wi;
}

static auto IncrementTargetCurr(WindowIndices wi) -> WindowIndices {
  ++wi.target_indices.curr;
  return wi;
}

static auto IncrementTargetCurr(WindowIndices wi, std::uint32_t const val)
    -> WindowIndices {
  wi.target_indices.curr += val;
  return wi;
}

static auto IncrementCurr(WindowIndices wi) -> WindowIndices {
  return IncrementTargetCurr(IncrementQueryCurr(wi));
}

static auto BindFirst(WindowIndices wi) -> WindowIndices {
  wi.query_indices.active_interval.first = wi.query_indices.curr;
  wi.target_indices.active_interval.first = wi.target_indices.curr;

  return wi;
}

static auto TryBindFirst(WindowIndices wi) -> WindowIndices {
  if (wi.query_indices.active_interval.first == kInvalidIndex) {
    return BindFirst(wi);
  }

  return wi;
}

static auto BindLast(WindowIndices wi) -> WindowIndices {
  wi.query_indices.active_interval.last = wi.query_indices.curr + 1;
  wi.target_indices.active_interval.last = wi.target_indices.curr + 1;

  return wi;
}

static auto TryBindLast(WindowIndices wi) -> WindowIndices {
  if (wi.query_indices.active_interval.first != kInvalidIndex) {
    return BindLast(wi);
  }

  return wi;
}

static auto GetWindowIdx(std::span<ReferenceWindow> windows,
                         std::uint32_t ovlp_start) -> std::uint32_t {
  return std::distance(
      windows.begin(),
      std::upper_bound(windows.begin(), windows.end(), ovlp_start,
                       [](std::uint32_t const ovlp_start,
                          ReferenceWindow const& ref_window) -> bool {
                         return ovlp_start < ref_window.interval.last;
                       }));
}

static auto MakeAlignedSegment(std::uint32_t const window_first,
                               Interval const target_interval,
                               Interval const query_interval,
                               NucleicView query_view,
                               double const quality_threshold)
    -> std::optional<AlignedSegment> {
  auto quality_sum = 0.0;
  for (auto pos = query_interval.first; pos < query_interval.last; ++pos) {
    quality_sum += query_view.Quality(pos) - 33;
  }

  if (auto const quality_avg = quality_sum / IntervalLength(query_interval);
      quality_avg < quality_threshold) {
    return std::nullopt;
  }

  return AlignedSegment{
      .alignment_local_interval =
          LocalizeInterval(window_first, target_interval),
      .bases = query_view.InflateData(query_interval.first,
                                      IntervalLength(query_interval)),
      .quality = query_view.InflateQuality(query_interval.first,
                                           IntervalLength(query_interval))};
}

static auto BindOverlapToWindow(NucleicView query_view,
                                biosoup::Overlap const& ovlp,
                                std::span<ReferenceWindow> windows,
                                std::uint32_t win_idx,
                                double const quality_threshold) -> void {
  auto cigar = std::string_view(ovlp.alignment);
  auto indices = MakeWindowIndices(ovlp, query_view.InflatedLenght());

  auto try_update_window = [&win_idx, query_view, windows, quality_threshold](
                               WindowIndices indices) -> WindowIndices {
    if (indices.target_indices.curr + 1 != windows[win_idx].interval.last) {
      return indices;
    }

    if (auto opt_segment =
            MakeAlignedSegment(windows[win_idx].interval.first,
                               indices.target_indices.active_interval,
                               indices.query_indices.active_interval,
                               query_view, quality_threshold);
        opt_segment) {
      windows[win_idx].aligned_segments.push_back(*opt_segment);
    };

    indices.query_indices.active_interval = {kInvalidIndex, kInvalidIndex};
    indices.target_indices.active_interval = {kInvalidIndex, kInvalidIndex};
    ++win_idx;

    return indices;
  };

  for (auto j = 0U, i = 0U; i < cigar.size(); ++i) {
    if (IsMisOrMatch(cigar[i])) {
      auto const n = std::atoi(&cigar[j]);
      j = i + 1;

      for (auto k = 0; k < n; ++k) {
        indices =
            try_update_window(BindLast(TryBindFirst(IncrementCurr(indices))));
      }
    } else if (IsInsertion(cigar[i])) {
      indices = IncrementQueryCurr(indices, std::atoi(&cigar[j]));
      j = i + 1;
    } else if (IsDeletion(cigar[i])) {
      auto const n = std::atoi(&cigar[j]);
      j = i + 1;

      for (auto k = 0; k < n; ++k) {
        indices = try_update_window(IncrementTargetCurr(indices));
      }
    } else if (IsClipOrPad(cigar[i])) {
      j = i + 1;
    }
  }
}

static auto BindReadSegmentsToWindows(
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> reads,
    std::span<biosoup::Overlap const> overlaps,
    std::span<ReferenceWindow> windows, double const quality_threshold)
    -> void {
  for (auto ovlp_idx = 0U; ovlp_idx < overlaps.size(); ++ovlp_idx) {
    auto win_idx = GetWindowIdx(windows, overlaps[ovlp_idx].rhs_begin);
    if (win_idx >= windows.size()) {
      continue;
    }

    BindOverlapToWindow(NucleicView(reads[overlaps[ovlp_idx].lhs_id].get(),
                                    !overlaps[ovlp_idx].strand),
                        overlaps[ovlp_idx], windows, win_idx,
                        quality_threshold);
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

auto MakePOAWindows(
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

auto MakeConsensus(std::string_view backbone_data,
                   std::string_view backbone_quality,
                   ReferenceWindowView ref_window_view, POAConfig poa_cfg)
    -> ConsensusResult {
  if (ref_window_view.aligned_segments.size() < 2) {
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
