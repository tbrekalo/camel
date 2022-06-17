#include "camel/correct.h"

#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <numeric>
#include <random>

#include "biosoup/timer.hpp"
#include "camel/io.h"
#include "camel/mapping.h"
#include "detail/overlap.h"
#include "edlib.h"
#include "fmt/core.h"
#include "fmt/ostream.h"
#include "fmt/ranges.h"
#include "spoa/alignment_engine.hpp"
#include "spoa/spoa.hpp"
#include "tsl/robin_map.h"

namespace camel {

namespace detail {

static auto constexpr kAllowedFuzzPercent = 0.01;
static auto constexpr kSmallWindowPercent = 0.04;

static auto constexpr kMinCoverage = 16U;

static auto constexpr kSnpProximityLimit = 20U;
static auto constexpr kCorrectBatchCap = 1UL << 24UL;

struct FastCovg {
  // mat, del, ins, mis
  std::array<std::uint_fast16_t, 4> signal;
};

struct Interval {
  std::uint32_t start_idx;
  std::uint32_t end_idx;
};

struct Evidence {
  std::uint32_t query_pos;
  char base;
};

// [[nodiscard]] static auto OverlapStrings(
//     tsl::robin_map<std::uint32_t, ReadOverlapsPair> const& reads_overlaps,
//     biosoup::Overlap const& ovlp) -> std::pair<std::string, std::string> {
//   auto query_str =
//       reads_overlaps.at(ovlp.lhs_id)
//           .read->InflateData(ovlp.lhs_begin, ovlp.lhs_end - ovlp.lhs_begin);
// 
//   auto target_str =
//       reads_overlaps.at(ovlp.rhs_id)
//           .read->InflateData(ovlp.rhs_begin, ovlp.rhs_end - ovlp.rhs_begin);
// 
//   if (!ovlp.strand) {
//     auto rc = biosoup::NucleicAcid("", target_str);
//     rc.ReverseAndComplement();
// 
//     target_str = rc.InflateData();
//   }
// 
//   return std::pair(std::move(query_str), std::move(target_str));
// }

[[nodiscard]] static auto AlignStrings(std::string const& query_str,
                                       std::string const& target_str)
    -> EdlibAlignResult {
  return edlibAlign(
      query_str.c_str(), query_str.size(), target_str.c_str(),
      target_str.size(),
      edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));
}

// [[nodiscard]] static auto FindAlignments(
//     State& state,
//     tsl::robin_map<std::uint32_t, ReadOverlapsPair> const& reads_overlaps,
//     std::uint32_t const query_id) -> std::vector<OverlapEdlibAlignment> {
//   auto dst = std::vector<OverlapEdlibAlignment>();
// 
//   auto const create_alignment =
//       [&reads_overlaps](biosoup::Overlap const& ovlp) -> OverlapEdlibAlignment {
//     auto const [query_str, target_str] = OverlapStrings(reads_overlaps, ovlp);
// 
//     return {.ovlp = ovlp, .edlib_result = AlignStrings(query_str, target_str)};
//   };
// 
//   dst.reserve(reads_overlaps.at(query_id).overlaps.size());
//   for (auto const& ovlp : reads_overlaps.at(query_id).overlaps) {
//     auto const rev_ovlp = ReverseOverlap(ovlp);
//     dst.push_back(create_alignment(rev_ovlp));
//   }
// 
//   for (auto const& [read_id, read_overlaps] : reads_overlaps) {
//     if (read_id != query_id) {
//       auto const& ovlps = read_overlaps.overlaps;
//       auto first = std::lower_bound(
//           ovlps.cbegin(), ovlps.cend(), query_id,
// 
//           [](biosoup::Overlap const& ovlp, std::uint32_t const lhs_id) -> bool {
//             return ovlp.lhs_id < lhs_id;
//           });
// 
//       for (; first != ovlps.cend() && first->lhs_id == query_id; ++first) {
//         dst.push_back(create_alignment(*first));
//       }
//     }
//   }
// 
//   return dst;
// }

// [[nodiscard]] static auto EstimateCoverage(
//     State& state,
//     tsl::robin_map<std::uint32_t, ReadOverlapsPair> const& reads_overlaps)
//     -> std::uint16_t {
//   auto const kSubsampleSize =
//       std::max(static_cast<std::size_t>(reads_overlaps.size() * 0.10),
//                std::min(100UL, reads_overlaps.size()));
// 
//   auto pile_futures = std::vector<std::future<std::uint16_t>>();
//   pile_futures.reserve(kSubsampleSize);
// 
//   // generate random ids
//   auto rng_ids = std::vector<std::uint32_t>(reads_overlaps.size());
//   std::transform(reads_overlaps.cbegin(), reads_overlaps.cend(),
//                  rng_ids.begin(),
//                  [](std::pair<std::uint32_t, ReadOverlapsPair> const& ro)
//                      -> std::uint32_t { return ro.first; });
// 
//   auto rng_engine = std::mt19937(42U);
//   std::shuffle(rng_ids.begin(), rng_ids.end(), rng_engine);
// 
//   rng_ids.resize(kSubsampleSize);
//   for (auto const rng_id : rng_ids) {
//     pile_futures.emplace_back(state.thread_pool->Submit(
//         [&reads_overlaps](std::uint32_t const read_id) -> std::uint16_t {
//           auto constexpr kShrinkShift = 3U;
// 
//           auto const& backbone_read = reads_overlaps.at(read_id).read;
//           auto pile = std::vector<std::uint16_t>(
//               ((backbone_read->inflated_len) >> kShrinkShift) + 2U, 0U);
// 
//           auto events = std::vector<std::uint32_t>();
//           for (auto const& ovlp : reads_overlaps.at(read_id).overlaps) {
//             events.push_back(((ovlp.rhs_begin >> kShrinkShift) + 1U) << 1U);
//             events.push_back((((ovlp.rhs_end >> kShrinkShift) - 1U) << 1U) |
//                              1U);
//           }
// 
//           for (auto const& [that_id, that_ro] : reads_overlaps) {
//             if (that_id != read_id) {
//               auto ovlp_iter =
//                   std::lower_bound(that_ro.overlaps.cbegin(),
//                                    that_ro.overlaps.cend(), backbone_read->id,
//                                    [](biosoup::Overlap const& ovlp,
//                                       std::uint32_t const read_id) -> bool {
//                                      return ovlp.lhs_id < read_id;
//                                    });
// 
//               for (; ovlp_iter != that_ro.overlaps.cend() &&
//                      ovlp_iter->lhs_id == backbone_read->id;
//                    ++ovlp_iter) {
//                 events.push_back(((ovlp_iter->lhs_begin >> kShrinkShift) + 1U)
//                                  << 1U);
//                 events.push_back(
//                     (((ovlp_iter->lhs_end >> kShrinkShift) - 1U) << 1U) | 1U);
//               }
//             }
//           }
// 
//           auto convg = 0U;
//           auto last_event = 0U;
//           std::sort(events.begin(), events.end());
//           for (auto const event : events) {
//             for (auto i = last_event; i < (event >> 1U); ++i) {
//               if (convg > 0U) {
//                 pile[i] =
//                     std::clamp(pile[i] + convg, 0U,
//                                1U * std::numeric_limits<std::uint16_t>::max());
//               }
//             }
// 
//             convg += 1 - (2 * (event & 1));
//             last_event = (event >> 1U);
//           }
// 
//           std::nth_element(pile.begin(),
//                            std::next(pile.begin(), pile.size() / 2U),
//                            pile.end());
// 
//           return pile[pile.size() / 2U];
//         },
//         rng_id));
//   }
// 
//   auto medians = std::vector<std::uint16_t>(pile_futures.size());
//   std::transform(pile_futures.begin(), pile_futures.end(), medians.begin(),
//                  std::mem_fn(&std::future<std::uint16_t>::get));
// 
//   std::nth_element(medians.begin(),
//                    std::next(medians.begin(), medians.size() / 2U),
//                    medians.end());
// 
//   return medians[medians.size() / 2U];
// }
// 

}  // namespace detail

auto SnpErrorCorrect(
    State& state, MapCfg const map_cfg, CorrectConfig const correct_cfg,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> src_reads)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
  auto dst = std::vector<std::unique_ptr<biosoup::NucleicAcid>>();

  auto timer = biosoup::Timer();

  timer.Start();
  // auto reads_overlaps = detail::FindOverlapsAndFilterReads(
  //     state, map_cfg, correct_cfg.correct_window, correct_cfg.depth,
  //     std::move(src_reads));
  // fmt::print(stderr,
  //            "[camel::SnpErrorCorrect]({:12.3f}) tied reads with overlaps\n",
  //            timer.Stop());

  // timer.Start();
  // auto const kEstimatedCovg = detail::EstimateCoverage(state, reads_overlaps);
  // fmt::print(stderr,
  //            "[camel::SnpErrorCorrect]({:12.3f}) estimated coverage: {}\n",
  //            timer.Stop(), kEstimatedCovg);


  fmt::print(stderr, "[camel::SnpErrorCorrect]({:12.3f})\n",
             timer.elapsed_time());

  return dst;
}

}  // namespace camel
