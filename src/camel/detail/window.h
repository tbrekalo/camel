#ifndef CAMEL_DETAIL_WINDOW_H_
#define CAMEL_DETAIL_WINDOW_H_

#include <memory>
#include <span>
#include <vector>

#include "alignment.h"
#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"
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

[[nodiscard]] auto CreateWindowsFromAlignments(
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> reads,
    std::span<biosoup::Overlap const> overlaps,
    std::span<EdlibAlignResult const> edlib_results,
    std::uint32_t const global_coverage_estimate)
    -> std::vector<ReferenceWindow>;

}  // namespace camel::detail

#endif /* CAMEL_DETAIL_WINDOW_H_ */
