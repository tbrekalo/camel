#include "alignment.h"

#include "nucleic_view.h"

namespace camel::detail {

auto ExtractSubstrings(
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> reads,
    biosoup::Overlap ovlp) -> std::tuple<std::string, std::string> {
  auto const lhs_view = NucleicView(reads[ovlp.lhs_id].get(), false);
  auto const rhs_view = NucleicView(reads[ovlp.rhs_id].get(),
                                    /* is_reverse_complement = */ !ovlp.strand);

  auto lhs_str =
      lhs_view.InflateData(ovlp.lhs_begin, ovlp.lhs_end - ovlp.lhs_begin);
  auto rhs_str =
      rhs_view.InflateData(ovlp.rhs_begin, ovlp.rhs_end - ovlp.rhs_begin);

  return {std::move(lhs_str), std::move(rhs_str)};
}

auto AlignStrings(std::string_view lhs_str_view, std::string_view rhs_str_view)
    -> EdlibAlignResult {
  /* clang-format off */
  return edlibAlign(
      lhs_str_view.data(), lhs_str_view.length(), 
      rhs_str_view.data(), rhs_str_view.length(),
      edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));
  /* clang-format on */
}

auto OverlapToALignment(
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> reads,
    biosoup::Overlap ovlp) -> EdlibAlignResult {
  return std::apply(AlignStrings, ExtractSubstrings(reads, ovlp));
}

}  // namespace camel::detail
