#ifndef CAMEL_DETAIL_ALIGNMENT_H_
#define CAMEL_DETAIL_ALIGNMENT_H_

#include <memory>
#include <span>

#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"

namespace camel::detail {

auto IsMisOrMatch(char chr) -> bool;
auto IsInsertion(char chr) -> bool;
auto IsDeletion(char chr) -> bool;
auto IsClipOrPad(char chr) -> bool;

[[nodiscard]] auto AlignedOverlap(
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> reads,
    biosoup::Overlap ovlp) -> biosoup::Overlap;

}  // namespace camel::detail

#endif /* CAMEL_DETAIL_ALIGNMENT_H_ */
