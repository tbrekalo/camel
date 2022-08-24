#ifndef CAMEL_DETAIL_COVERAGE_H_
#define CAMEL_DETAIL_COVERAGE_H_

#include <algorithm>
#include <array>
#include <cstdint>
#include <memory>
#include <vector>

#include "alignment.h"
#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"
#include "nonstd/span.hpp"
#include "tbb/task_arena.h"

namespace camel::detail {

struct CoverageSignals {
  std::array<std::uint16_t, 6U> signals;

  static inline auto constexpr kDelIdx = 4U;
  static inline auto constexpr kInsIdx = 5U;
};

[[nodiscard]] auto EstimateCoverage(
    tbb::task_arena& task_arena,
    nonstd::span<std::unique_ptr<biosoup::NucleicAcid>> reads,
    nonstd::span<std::vector<biosoup::Overlap>> overlaps) -> std::uint16_t;

[[nodiscard]] auto CalculateCoverage(
    nonstd::span<std::unique_ptr<biosoup::NucleicAcid>> reads,
    nonstd::span<biosoup::Overlap> overlaps,
    nonstd::span<EdlibAlignResult> edlib_results)
    -> std::vector<CoverageSignals>;
}  // namespace camel::detail

#endif /* CAMEL_DETAIL_COVERAGE_H_ */
