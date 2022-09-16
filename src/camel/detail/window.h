#ifndef CAMEL_DETAIL_WINDOW_H_
#define CAMEL_DETAIL_WINDOW_H_

#include <memory>
#include <span>
#include <vector>

#include "alignment.h"
#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"
#include "camel/poa_config.h"
#include "interval.h"

namespace camel::detail {

inline auto constexpr kWinPadding = 13U;
inline auto constexpr kAllowedFuzzPercent = 0.05;
inline auto constexpr kSmallWindowPercent = 0.05;

struct AlignedSegment {
  Interval alignment_local_interval;
  std::string bases;
};

struct ReferenceWindow {
  Interval interval;
  std::vector<AlignedSegment> aligned_segments;
};

struct ReferenceWindowView {
  Interval interval;
  std::span<AlignedSegment> aligned_segments;
};

[[nodiscard]] auto CreateWindowsFromAlignments(
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> reads,
    std::span<biosoup::Overlap const> overlaps,
    std::span<EdlibAlignResult const> edlib_results,
    std::uint32_t const target_window_len,
    std::uint32_t const global_coverage_estimate)
    -> std::vector<ReferenceWindow>;

auto ReleaseAlignmentEngines() -> std::size_t;

[[nodiscard]] auto WindowConsensus(std::string_view backbone_view,
                                   ReferenceWindowView, POAConfig poa_cfg)
    -> std::string;

}  // namespace camel::detail

#endif /* CAMEL_DETAIL_WINDOW_H_ */
