#include "alignment.h"

#include "bindings/cpp/WFAligner.hpp"
#include "edlib.h"
#include "nucleic_view.h"
#include "tbb/enumerable_thread_specific.h"

namespace camel::detail {

auto IsMisOrMatch(char chr) -> bool {
  return chr == 'M' || chr == '=' || chr == 'X';
}

auto IsInsertion(char chr) -> bool { return chr == 'I'; }

auto IsDeletion(char chr) -> bool { return chr == 'D' || chr == 'N'; }

auto IsClipOrPad(char chr) -> bool {
  return chr == 'S' || chr == 'H' || chr == 'P';
}

static auto Aligners =
    tbb::enumerable_thread_specific<std::unique_ptr<wfa::WFAlignerGapAffine>>();

static auto GetAligner(std::int8_t mismatch, std::int8_t grap_opening,
                       std::int8_t gap_extension)
    -> std::unique_ptr<wfa::WFAlignerGapAffine>& {
  auto& aligner = Aligners.local();
  if (!aligner) {
    aligner = std::make_unique<wfa::WFAlignerGapAffine>(
        mismatch, grap_opening, gap_extension, wfa::WFAligner::Alignment,
        wfa::WFAligner::MemoryUltralow);

    aligner->setHeuristicWFadaptive(10, 50, 1);
    aligner->setHeuristicBandedAdaptive(-50, +50, 1);
    aligner->setHeuristicXDrop(100, 100);
    aligner->setHeuristicZDrop(100, 100);
  }
  return aligner;
};

[[nodiscard]] auto ExtractSubstrings(
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

auto AlignedOverlap(
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> reads,
    biosoup::Overlap ovlp) -> biosoup::Overlap {
  auto [query_str, target_str] = ExtractSubstrings(reads, ovlp);
  auto& aligner = GetAligner(5, 4, 2);
  aligner->alignEnd2End(query_str, target_str);
  ovlp.alignment = aligner->getAlignmentCigar();

  return ovlp;
}

}  // namespace camel::detail
