#include "alignment.h"

#include "nucleic_view.h"

namespace camel::detail {

auto ExtractSubstrings(
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> reads,
    biosoup::Overlap ovlp) -> std::tuple<std::string, std::string> {
  auto const query_view =
      NucleicView(reads[ovlp.lhs_id].get(),
                  /* is_reverse_complement = */ !ovlp.strand);
  auto const target_view = NucleicView(reads[ovlp.rhs_id].get(), false);

  auto query_str = query_view.InflateData(
      ovlp.strand ? ovlp.lhs_begin
                  : reads[ovlp.lhs_id]->inflated_len - ovlp.lhs_end,
      ovlp.lhs_end - ovlp.lhs_begin);

  auto target_str =
      target_view.InflateData(ovlp.rhs_begin, ovlp.rhs_end - ovlp.rhs_begin);

  return {std::move(query_str), std::move(target_str)};
}

auto AlignStrings(std::string_view query_str_view,
                  std::string_view target_str_view) -> EdlibAlignResult {
  /* clang-format off */
  return edlibAlign(
      query_str_view.data(), query_str_view.length(), 
      target_str_view.data(), target_str_view.length(),
      edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));
  /* clang-format on */
}

auto OverlapToALignment(
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> reads,
    biosoup::Overlap ovlp) -> EdlibAlignResult {
  return std::apply(AlignStrings, ExtractSubstrings(reads, ovlp));
}

}  // namespace camel::detail
