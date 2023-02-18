#include "coverage.h"

#include "tbb/parallel_for.h"

namespace camel::detail {

static auto constexpr kShrinkShift = 3U;

[[nodiscard]] static auto ReadMedianCoverageEstimate(
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> reads,
    std::span<std::vector<biosoup::Overlap> const> overlaps,
    std::uint32_t const read_id) -> std::uint16_t {
  auto pile = std::vector<std::uint16_t>(
      (reads[read_id]->inflated_len >> kShrinkShift) + 2U);

  auto events = std::vector<std::uint32_t>();
  events.reserve(2 * overlaps[read_id].size());

  for (auto const& ovlp : overlaps[read_id]) {
    events.push_back(((ovlp.rhs_begin >> kShrinkShift) + 1U) << 1U);
    events.push_back((((ovlp.rhs_end >> kShrinkShift) - 1U) << 1U) | 1U);
  }

  auto covg = 0U;
  auto last_event = 0U;
  std::sort(events.begin(), events.end());
  for (auto const event : events) {
    for (auto i = last_event; i < (event >> 1U); ++i) {
      pile[i] = std::clamp(pile[i] + covg, 0U,
                           1U * std::numeric_limits<std::uint16_t>::max());
    }

    covg += 1 - (2 * (event & 1));
    last_event = event >> 1U;
  }

  std::nth_element(pile.begin(), pile.begin() + pile.size() / 2, pile.end());
  return pile[pile.size() / 2];
}

auto EstimateCoverage(
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> reads,
    std::span<std::vector<biosoup::Overlap> const> overlaps) -> std::uint16_t {
  auto covg_medians = std::vector<std::uint16_t>(reads.size());
  tbb::parallel_for(0UL, reads.size(), [&](std::size_t read_idx) -> void {
    covg_medians[read_idx] =
        detail::ReadMedianCoverageEstimate(reads, overlaps, read_idx);
  });

  std::nth_element(covg_medians.begin(),
                   covg_medians.begin() + covg_medians.size() / 2,
                   covg_medians.end());

  return covg_medians[covg_medians.size() / 2];
}

auto CalculateCoverage(
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> reads,
    std::span<biosoup::Overlap const> overlaps)
    -> std::vector<CoverageSignals> {
  auto const target_id = overlaps.front().rhs_id;
  auto dst = std::vector<CoverageSignals>(reads[target_id]->inflated_len + 1U);

  // for (auto i = 0U; i < reads[target_id]->inflated_len; ++i) {
  //   ++dst[i].val[reads[target_id]->Code(i)];
  // }

  // for (auto ovlp_idx = 0U; ovlp_idx < overlaps.size(); ++ovlp_idx) {
  //   auto const& edlib_res = edlib_results[ovlp_idx];
  //   auto const query_view =
  //       NucleicView(reads[overlaps[ovlp_idx].lhs_id].get(),
  //                   /* is_reverse_complement = */ !overlaps[ovlp_idx].strand);

  //   auto query_pos = overlaps[ovlp_idx].strand
  //                        ? overlaps[ovlp_idx].lhs_begin
  //                        : reads[overlaps[ovlp_idx].lhs_id]->inflated_len -
  //                              overlaps[ovlp_idx].lhs_end;
  //   auto target_pos = overlaps[ovlp_idx].rhs_begin;

  //   auto i = 0;
  //   for (; i < edlib_res.alignmentLength; ++i) {
  //     switch (edlib_res.alignment[i]) {
  //       case 0:  // match
  //       case 3:  // mismatch
  //         ++dst[target_pos].val[query_view.Code(query_pos)];
  //         ++query_pos;
  //         ++target_pos;
  //         break;
  //       case 1:  // insertion on the target
  //         ++dst[target_pos].val[CoverageSignals::kInsIdx];
  //         ++query_pos;
  //         break;
  //       case 2:  // deletion on the target
  //         ++dst[target_pos].val[CoverageSignals::kDelIdx];
  //         ++target_pos;
  //         break;
  //       default:
  //         break;
  //     }
  //   }
  // }

  return dst;
}

}  // namespace camel::detail
