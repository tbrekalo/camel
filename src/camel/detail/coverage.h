#ifndef CAMEL_DETAIL_COVERAGE_H_
#define CAMEL_DETAIL_COVERAGE_H_

#include <algorithm>
#include <array>
#include <cstdint>
#include <memory>
#include <vector>

#include "alignment.h"
#include "biosoup/overlap.hpp"
#include "nucleic_view.h"
#include "tbb/task_arena.h"

namespace camel::detail {

struct CoverageSignals {
  std::array<std::uint16_t, 6U> val;

  static inline auto constexpr kDelIdx = 4U;
  static inline auto constexpr kInsIdx = 5U;
};

[[nodiscard]] auto EstimateCoverage(
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> reads,
    std::span<std::vector<biosoup::Overlap> const> overlaps) -> std::uint16_t;

[[nodiscard]] auto CalculateCoverage(
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> reads,
    std::span<biosoup::Overlap const> overlaps,
    std::span<EdlibAlignResult const> edlib_results)
    -> std::vector<CoverageSignals>;
}  // namespace camel::detail

#endif /* CAMEL_DETAIL_COVERAGE_H_ */
