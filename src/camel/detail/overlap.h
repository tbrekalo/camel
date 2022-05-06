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
  kInvalid,
  kUnclassified,

  kInternal,      ///< overlap has a deep location in one of the sequnces
  kLhsContained,  ///< rhs sequence spans over lhs sequence
  kRhsContained,  ///< lhs sequence spans over rhs sequence
  kLhsToRhs,      ///< lhs sequence connects to rhs sequence
  kRhsToLhs,      ///< rhs sequence connects to lhs sequence

};

/**
 * @brief concept for overlap filter callable:
 *         Filter(biosoup::Overlap const&) -> bool;
 */
template <class T>
using IsOverlapFilter =
  std::is_invocable_r<bool, T, biosoup::Overlap const&>;

/**
 * @brief helper for @ref IsOverlapFilter
 */
template <class T>
inline bool constexpr IsOverlapFilerV = IsOverlapFilter<T>::value;

/**
 * @brief concept for overlap compare callable:
 *         OvlpCmp(biosoup::Overlap const&, biosoup::Overlap const&) -> bool 
 */
template <class T>
using IsOverlapCmpCallable =
    std::is_invocable_r<bool, T, biosoup::Overlap const&,
                        biosoup::Overlap const&>;

/**
 * @brief helper for @ref IsOverlapCmpCallable
 */
template <class T>
inline bool constexpr IsOverlapCmpCallableV =
  IsOverlapCmpCallable<T>::value;

/**
 * @brief concept for overlap sink callable:
 *          Sink(biosoup::Overlap const&) -> void;
 */
template <class T>
using IsOverlapSink = std::is_invocable_r<void, T, biosoup::Overlap const&>;

/**
 * @brief helper for @ref IsOverlapSink
 */
template <class T>
inline bool constexpr IsOverlapSinkV = IsOverlapSink<T>::value;

/**
 * @brief Determine overlap type provided origin sequence length information
 */
auto DetermineOverlapType(biosoup::Overlap const ovlp,
                          std::uint32_t const lhs_seq_size,
                          std::uint32_t const rhs_seq_size) -> OverlapType;

/**
 * @brief swap rhs and lhs
 */
auto ReverseOverlap(biosoup::Overlap const& ovlp) -> biosoup::Overlap;

/**
 * @brief calculate overlap length
 */
auto OverlapLength(biosoup::Overlap const& ovlp) -> std::uint32_t;

/**
 * @brief calculate overlap score
 */

auto OverlapScore(biosoup::Overlap const& ovlp) -> double;

}  // namespace camel::detail

#endif /* CAMEL_DETAIL_OVERLAP_H_ */
