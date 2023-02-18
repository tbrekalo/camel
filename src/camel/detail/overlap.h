#ifndef CAMEL_DETAIL_OVERLAP_H_
#define CAMEL_DETAIL_OVERLAP_H_

#include <cstdint>
#include <type_traits>

#include "biosoup/overlap.hpp"

namespace camel::detail {

/**
 * @brief specifies overlap types
 */
enum class OverlapType : std::uint8_t {
  kUnclassified,
  kInternal,      ///< overlap has a deep location in one of the sequnces
  kLhsContained,  ///< rhs sequence spans over lhs sequence
  kRhsContained,  ///< lhs sequence spans over rhs sequence
  kLhsToRhs,      ///< lhs sequence connects to rhs sequence
  kRhsToLhs,      ///< rhs sequence connects to lhs sequence

};

/**
 * @brief Determine overlap type provided origin sequence length information
 */
auto DetermineOverlapType(biosoup::Overlap const& ovlp,
                          std::uint32_t const lhs_seq_size,
                          std::uint32_t const rhs_seq_size) -> OverlapType;

auto OverlapLength(biosoup::Overlap const& ovlp) -> std::uint32_t;
auto OverlapError(biosoup::Overlap const& ovlp) -> double;
auto OverlapScore(biosoup::Overlap const& ovlp) -> double;

}  // namespace camel::detail

#endif /* CAMEL_DETAIL_OVERLAP_H_ */
