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
    aligner->setHeuristicBandedAdaptive(-10, +10, 1);
    aligner->setHeuristicXDrop(100, 100);
    aligner->setHeuristicZDrop(100, 100);
  }
  return aligner;
};

static auto ExtractSubstrings(
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

static auto WFAling(std::string& query, std::string& target) -> std::string {
  auto& aligner = GetAligner(5, 4, 2);
  aligner->alignEnd2End(query, target);
  return aligner->getAlignmentCigar();
}

static auto EdlibAlign(std::string_view lhs, std::string_view rhs)
    -> std::string {
  auto dst = std::string{};
  auto edlibRes = edlibAlign(
      lhs.data(), lhs.length(), rhs.data(), rhs.length(),
      edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));

  if (edlibRes.status == EDLIB_STATUS_OK) {
    auto cigar_ptr = std::unique_ptr<char, void (*)(char*)>(
        edlibAlignmentToCigar(edlibRes.alignment, edlibRes.alignmentLength,
                              EDLIB_CIGAR_STANDARD),
        +[](char* chr) -> void { free(chr); });
    dst = std::string(cigar_ptr.get());
  }

  edlibFreeAlignResult(edlibRes);
  return dst;
}

auto AlignedOverlap(
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> reads,
    biosoup::Overlap ovlp) -> biosoup::Overlap {
  auto [query_str, target_str] = ExtractSubstrings(reads, ovlp);
  ovlp.alignment = EdlibAlign(query_str, target_str);

  return ovlp;
}

}  // namespace camel::detail
