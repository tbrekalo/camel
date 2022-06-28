#include "camel/correct.h"

#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <numeric>
#include <random>
#include <variant>

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

static auto constexpr kShrinkShift = 3U;
static auto constexpr kMaxInsLen = 32U;

struct MatchInterval {
  std::uint32_t first;
  std::uint32_t last;
};

struct MismatchSite {
  std::uint32_t pos;
};

struct DeletinInterval {
  std::uint32_t first;
  std::uint32_t last;
};

struct InsertionSite {
  std::uint32_t pos;
  std::uint64_t val;
  std::uint8_t len;
};

using AlignmentUnit =
    std::variant<MatchInterval, DeletinInterval, InsertionSite>;

class EdlibAlignRAIIRes : public EdlibAlignResult {
  bool engaged = false;

 public:
  auto Release() noexcept {
    if (engaged) {
      edlibFreeAlignResult(*static_cast<EdlibAlignResult*>(this));
      engaged = false;
    }
  }

  EdlibAlignRAIIRes(EdlibAlignResult const& eres) noexcept
      : EdlibAlignResult(eres), engaged(true) {}

  EdlibAlignRAIIRes& operator=(EdlibAlignResult const& eres) noexcept {
    Release();

    *this = EdlibAlignRAIIRes(eres);
    return *this;
  }

  EdlibAlignRAIIRes(EdlibAlignRAIIRes const&) = delete;
  EdlibAlignRAIIRes& operator=(EdlibAlignRAIIRes const&) = delete;

  EdlibAlignRAIIRes(EdlibAlignRAIIRes&& that) noexcept {
    *this = std::move(that);
  }

  EdlibAlignRAIIRes& operator=(EdlibAlignRAIIRes&& that) noexcept {
    Release();

    engaged = std::exchange(that.engaged, false);

    status = that.status;
    editDistance = that.editDistance;
    endLocations = that.endLocations;
    startLocations = that.startLocations;
    numLocations = that.numLocations;
    alignment = that.alignment;
    alignmentLength = that.alignmentLength;
    alphabetLength = that.alphabetLength;

    return *this;
  }

  ~EdlibAlignRAIIRes() noexcept { Release(); }
};

[[nodiscard]] static auto EstimateCoverage(
    State& state,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads,
    std::vector<std::vector<biosoup::Overlap>> const& overlaps)
    -> std::uint32_t {
  auto covg_futures = std::vector<std::future<std::uint16_t>>();
  covg_futures.reserve(reads.size());

  auto const covg_task =
      [&reads, &overlaps](std::uint32_t const read_id) -> std::uint16_t {
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
  };

  for (auto read_id = 0U; read_id < reads.size(); ++read_id) {
    if (!overlaps[read_id].empty()) {
      covg_futures.emplace_back(state.thread_pool->Submit(covg_task, read_id));
    }
  }

  auto covgs = std::vector<std::uint16_t>(covg_futures.size());
  std::transform(covg_futures.begin(), covg_futures.end(), covgs.begin(),
                 std::mem_fn(&std::future<std::uint16_t>::get));

  std::nth_element(covgs.begin(), covgs.begin() + covgs.size() / 2,
                   covgs.end());

  return covgs[covgs.size() / 2];
}

[[nodiscard]] static auto ExtractSubstrings(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads,
    biosoup::Overlap const& ovlp) -> std::pair<std::string, std::string> {
  auto lhs_str = reads[ovlp.lhs_id]->InflateData(ovlp.lhs_begin,
                                                 ovlp.lhs_end - ovlp.rhs_begin);
  auto rhs_str = reads[ovlp.rhs_id]->InflateData(ovlp.rhs_begin,
                                                 ovlp.rhs_end - ovlp.rhs_begin);

  if (!ovlp.strand) {
    auto acid = biosoup::NucleicAcid("", rhs_str);
    acid.ReverseAndComplement();

    rhs_str = acid.InflateData();
  }

  return {lhs_str, rhs_str};
};

[[nodiscard]] static auto AlignStrings(std::string_view lhs_str_view,
                                       std::string_view rhs_str_view)
    -> EdlibAlignRAIIRes {
  return edlibAlign(lhs_str_view.data(), lhs_str_view.length(),
                    rhs_str_view.data(), rhs_str_view.length(),
                    edlibDefaultAlignConfig());
}

[[nodiscard]] static auto CompressedAlignment(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads,
    biosoup::Overlap const& ovlp) -> std::vector<AlignmentUnit> {
  auto [lhs_substr, rhs_substr] = ExtractSubstrings(reads, ovlp);
  auto const edlib_res = AlignStrings(lhs_substr, rhs_substr);
  auto dst = std::vector<AlignmentUnit>();

  dst.reserve(edlib_res.alignmentLength / 10U);

  auto lhs_pos = ovlp.lhs_begin;
  auto rhs_pos = 0U;
  for (auto i = 0U; i < edlib_res.alignmentLength;) {
    switch (edlib_res.alignment[i]) {
      case 0: {
        auto start = lhs_pos;
        while (i < edlib_res.alignmentLength && edlib_res.alignment[i] == 0) {
          ++lhs_pos;
          ++rhs_pos;
          ++i;
        }

        dst.emplace_back(MatchInterval{.first = start, .last = lhs_pos});
        break;
      }
      case 1: {
        auto start = lhs_pos;
        while (i < edlib_res.alignmentLength && edlib_res.alignment[i] == 1) {
          ++lhs_pos;
          ++i;
        }

        dst.emplace_back(DeletinInterval{.first = start, .last = lhs_pos});
        break;
      }

      case 2: {
        auto storage = 0ULL;
        auto len = 0U;
        while (i < edlib_res.alignmentLength && edlib_res.alignment[i] == 2) {
          storage = len < kMaxInsLen
                        ? ((storage << 2 |
                            biosoup::kNucleotideCoder[rhs_substr[rhs_pos]]))
                        : storage;
          len += (len < kMaxInsLen);

          ++rhs_pos;
          ++i;
        }

        dst.emplace_back(InsertionSite{.pos = lhs_pos,
                                       .val = storage,
                                       .len = static_cast<std::uint8_t>(len)});
        break;
      }
      case 3: {
        auto start = lhs_pos;
        while (i < edlib_res.alignmentLength && edlib_res.alignment[i] == 0) {
          ++lhs_pos;
          ++rhs_pos;
          ++i;
        }

        dst.emplace_back(DeletinInterval{.first = start, .last = lhs_pos});
        break;
      }

      default:
        ++i;
        break;
    }
  }

  return dst;
}

[[nodiscard]] static auto GenerateConsensus(
    std::unique_ptr<biosoup::NucleicAcid> const& read,
    std::vector<std::vector<AlignmentUnit>> const& alignments)
    -> std::unique_ptr<biosoup::NucleicAcid> {
  return nullptr;
}

}  // namespace detail

auto ErrorCorrect(State& state, MapCfg const map_cfg,
                  CorrectConfig const correct_cfg,
                  std::vector<std::unique_ptr<biosoup::NucleicAcid>> src_reads)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
  auto dst = std::vector<std::unique_ptr<biosoup::NucleicAcid>>();

  auto timer = biosoup::Timer();

  timer.Start();
  auto const ovlps = FindOverlaps(state, map_cfg, src_reads);

  auto target_ids = std::vector<std::uint32_t>();
  for (auto read_id = 0U; read_id < ovlps.size(); ++read_id) {
    if (!ovlps[read_id].empty()) {
      target_ids.push_back(read_id);
    }
  }

  auto const kNTargets = target_ids.size();
  auto const kNContained = ovlps.size() - kNTargets;

  timer.Stop();

  timer.Start();
  auto const kCovgEstimate = detail::EstimateCoverage(state, src_reads, ovlps);
  fmt::print(stderr, "[camel::ErrorCorrect]({:12.3f}) coverage estimate: {}\n",
             timer.Stop(), kCovgEstimate);

  {
    timer.Start();
    auto alignments =
        std::vector<std::vector<std::vector<detail::AlignmentUnit>>>(kNTargets);

    auto align_futures =
        std::vector<std::future<std::vector<detail::AlignmentUnit>>>();
    align_futures.reserve(kNTargets * 32U);

    for (auto i = 0U; i < kNTargets; ++i) {
      auto const target_id = target_ids[i];
      if (!ovlps[target_id].empty()) {
        alignments[i].reserve(ovlps[target_id].size());
        for (auto const& ovlp : ovlps[target_id]) {
          align_futures.emplace_back(
              state.thread_pool->Submit(detail::CompressedAlignment,
                                        std::cref(src_reads), std::cref(ovlp)));
        }
      }
    }

    auto consensus_futures =
        std::vector<std::future<std::unique_ptr<biosoup::NucleicAcid>>>();
    consensus_futures.reserve(kNTargets);

    auto align_future_id = 0U;
    for (auto i = 0U; i < kNTargets; ++i) {
      auto const target_id = target_ids[i];
      for (auto ovlp_id = 0U; ovlp_id < ovlps[target_id].size(); ++ovlp_id) {
        alignments[i].push_back(align_futures[align_future_id++].get());
      }

      consensus_futures.emplace_back(state.thread_pool->Submit(
          [](std::unique_ptr<biosoup::NucleicAcid>& read,
             std::vector<std::vector<detail::AlignmentUnit>>& alignments)
              -> std::unique_ptr<biosoup::NucleicAcid> {
            auto dst = detail::GenerateConsensus(read, alignments);

            read.reset();
            alignments.clear();

            return dst;
          },
          std::ref(src_reads[target_id]), std::ref(alignments[i])));

      if ((i & 127U) == 0 || i + 1U == kNTargets) {
        fmt::print(stderr,
                   "\r[camel::ErrorCorrect]({:12.3f}) collected {} / {} "
                   "alignment futures",
                   timer.Lap(), i, kNTargets);
      }
    }

    timer.Stop();
    fmt::print(stderr, "\n");
  }

  fmt::print(stderr, "[camel::ErrorCorrect]({:12.3f})\n", timer.elapsed_time());

  return dst;
}

}  // namespace camel
