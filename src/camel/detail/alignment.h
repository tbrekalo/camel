#ifndef CAMEL_DETAIL_ALIGNMENT_H_
#define CAMEL_DETAIL_ALIGNMENT_H_

#include <memory>

#include "biosoup/overlap.hpp"
#include "edlib.h"
#include "nucleic_view.h"

namespace camel::detail {

[[nodiscard]] auto ExtractSubstrings(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads,
    biosoup::Overlap const& ovlp) -> std::pair<std::string, std::string>;

[[nodiscard]] static auto AlignStrings(std::string_view lhs_str_view,
                                       std::string_view rhs_str_view)
    -> EdlibAlignResult;

}  // namespace camel::detail

#endif /* CAMEL_DETAIL_ALIGNMENT_H_ */
