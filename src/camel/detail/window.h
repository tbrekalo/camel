#ifndef CAMEL_DETAIL_WINDOW_H_
#define CAMEL_DETAIL_WINDOW_H_

#include <memory>
#include <span>
#include <vector>

#include "alignment.h"
#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"
#include "camel/correct.h"
#include "interval.h"

namespace camel::detail {

inline auto constexpr kWinPadding = 13U;
inline auto constexpr kAllowedFuzzPercent = 0.01;
inline auto constexpr kSmallWindowPercent = 0.05;

struct AlignedSegment {
  Interval alignment_local_interval;
  std::string bases;
  std::string quality;
};

struct ReferenceWindow {
  Interval interval;
  std::vector<AlignedSegment> aligned_segments;
};

struct ReferenceWindowView {
  Interval interval;
  std::span<AlignedSegment const> aligned_segments;
};

struct ConsensusResult {
  std::string bases;
  bool is_corrected;
};

[[nodiscard]] auto CreateWindowsFromAlignments(
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> reads,
    std::span<biosoup::Overlap const> overlaps,
    WindowConfig window_cfg,
    std::uint32_t const global_coverage_estimate)
    -> std::vector<ReferenceWindow>;

auto ReleaseAlignmentEngines() -> std::size_t;

[[nodiscard]] auto WindowConsensus(std::string_view backbone_data,
                                   std::string_view backbone_quality,
                                   ReferenceWindowView, POAConfig poa_cfg)
    -> ConsensusResult;

}  // namespace camel::detail

#endif /* CAMEL_DETAIL_WINDOW_H_ */
