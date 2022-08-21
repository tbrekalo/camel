#ifndef CAMEL_DETAIL_COVERAGE_H_
#define CAMEL_DETAIL_COVERAGE_H_

#include <algorithm>
#include <array>
#include <cstdint>
#include <memory>
#include <vector>

#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"
#include "alignment.h"

namespace camel::detail {

struct CoverageSignals {
  std::array<std::uint16_t, 6U> signals;

  static inline auto constexpr kDelIdx = 4U;
  static inline auto constexpr kInsIdx = 5U;
};

[[nodiscard]] auto ReadMedianCoverageEstimate(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads,
    std::vector<std::vector<biosoup::Overlap>> const& overlaps,
    std::uint32_t const read_id) -> std::uint16_t;

[[nodiscard]] auto CalculateCoverage(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads,
    std::vector<biosoup::Overlap> const& overlaps,
    std::vector<EdlibAlignResult> const& edlib_results)
    -> std::vector<CoverageSignals>;
}  // namespace camel::detail

#endif /* CAMEL_DETAIL_COVERAGE_H_ */
