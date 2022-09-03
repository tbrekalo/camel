#ifndef CAMEL_DETAIL_ALIGNMENT_H_
#define CAMEL_DETAIL_ALIGNMENT_H_

#include <memory>

#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"
#include "edlib.h"
#include "nonstd/span.hpp"

namespace camel::detail {

[[nodiscard]] auto ExtractSubstrings(
    nonstd::span<std::unique_ptr<biosoup::NucleicAcid>> reads,
    biosoup::Overlap ovlp) -> std::tuple<std::string, std::string>;

[[nodiscard]] auto AlignStrings(std::string_view lhs_str_view,
                                std::string_view rhs_str_view)
    -> EdlibAlignResult;

[[nodiscard]] auto OverlapToALignment(
    nonstd::span<std::unique_ptr<biosoup::NucleicAcid>> reads,
    biosoup::Overlap ovlp) -> EdlibAlignResult;

}  // namespace camel::detail

#endif /* CAMEL_DETAIL_ALIGNMENT_H_ */
