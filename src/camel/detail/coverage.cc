#include "coverage.h"

#include "tbb/parallel_for.h"

namespace camel::detail {

static auto constexpr kShrinkShift = 3U;

[[nodiscard]] static auto ReadMedianCoverageEstimate(
    nonstd::span<std::unique_ptr<biosoup::NucleicAcid>> const& reads,
    nonstd::span<std::vector<biosoup::Overlap>> const& overlaps,
    std::uint32_t const read_id) -> std::uint16_t {
  auto pile = std::vector<std::uint16_t>(
      (reads[read_id]->inflated_len >> kShrinkShift) + 2U);

  auto events = std::vector<std::uint32_t>();
  events.reserve(2 * overlaps[read_id].size());

  for (auto const& ovlp : overlaps[read_id]) {
    events.push_back(((ovlp.lhs_begin >> kShrinkShift) + 1U) << 1U);
    events.push_back((((ovlp.lhs_end >> kShrinkShift) - 1U) << 1U) | 1U);
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

auto EstimateCoverage(nonstd::span<std::unique_ptr<biosoup::NucleicAcid>> reads,
                      nonstd::span<std::vector<biosoup::Overlap>> overlaps)
    -> std::uint16_t {
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
    nonstd::span<std::unique_ptr<biosoup::NucleicAcid>> reads,
    nonstd::span<biosoup::Overlap> overlaps,
    nonstd::span<EdlibAlignResult> edlib_results)
    -> std::vector<CoverageSignals> {
  auto const query_id = overlaps.front().lhs_id;
  auto dst = std::vector<CoverageSignals>(reads[query_id]->inflated_len + 1U);

  for (auto i = 0U; i < reads[query_id]->inflated_len; ++i) {
    ++dst[i].signals[reads[query_id]->Code(i)];
  }

  for (auto ovlp_idx = 0U; ovlp_idx < overlaps.size(); ++ovlp_idx) {
    auto const& edlib_res = edlib_results[ovlp_idx];
    auto const rhs_view = NucleicView(reads[overlaps[ovlp_idx].rhs_id].get(),
                                      !overlaps[ovlp_idx].strand);

    auto lhs_pos = overlaps[ovlp_idx].lhs_begin;
    auto rhs_pos = overlaps[ovlp_idx].rhs_begin;

    auto i = 0;
    for (; i < edlib_res.alignmentLength && edlib_res.alignment[i] == 2; ++i)
      ;  // skip initial insertion
    for (; i < edlib_res.alignmentLength; ++i) {
      switch (edlib_res.alignment[i]) {
        case 0:  // match
        case 3:  // mismatch
          ++dst[lhs_pos].signals[rhs_view.Code(rhs_pos)];
          ++lhs_pos;
          ++rhs_pos;
          break;
        case 1:  // deletion on the query
          ++dst[lhs_pos].signals[CoverageSignals::kDelIdx];
          ++lhs_pos;
          break;
        case 2:  // insertion on the query
          ++dst[lhs_pos - 1U].signals[CoverageSignals::kInsIdx];
          ++rhs_pos;
          break;
        default:
          break;
      }
    }
  }

  return dst;
}

}  // namespace camel::detail
