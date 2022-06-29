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
#include "ram/minimizer_engine.hpp"
#include "spoa/alignment_engine.hpp"
#include "spoa/spoa.hpp"
#include "tsl/robin_map.h"

namespace camel {

namespace detail {

static auto constexpr kShrinkShift = 3U;
static auto constexpr kMaxInsLen = 32U;

struct Insertion {
  std::uint64_t val;
  std::uint8_t len;
};

struct MatchInterval {
  std::uint32_t first;
  std::uint32_t last;
};

struct MismatchSite {
  std::uint32_t pos;
  std::uint8_t code;
};

struct DeletinInterval {
  std::uint32_t first;
  std::uint32_t last;
};

struct InsertionSite {
  std::uint32_t pos;
  Insertion ins;
};

struct AlignmentSummary {
  std::vector<MatchInterval> mats;
  std::vector<DeletinInterval> dels;
  std::vector<InsertionSite> inss;
  std::vector<MismatchSite> miss;
};

struct CoverageSignals {
  std::array<std::uint16_t, 6U> signals;

  static constexpr auto kDelIdx = 4U;
  static constexpr auto kInsIdx = 5U;
};

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
                                                 ovlp.lhs_end - ovlp.lhs_begin);
  auto rhs_str = reads[ovlp.rhs_id]->InflateData(ovlp.rhs_begin,
                                                 ovlp.rhs_end - ovlp.rhs_begin);

  if (!ovlp.strand) {
    auto acid = biosoup::NucleicAcid("", rhs_str);
    acid.ReverseAndComplement();

    rhs_str = acid.InflateData();
  }

  return {std::move(lhs_str), std::move(rhs_str)};
}

[[nodiscard]] static auto AlignStrings(std::string_view lhs_str_view,
                                       std::string_view rhs_str_view)
    -> EdlibAlignRAIIRes {
  return edlibAlign(
      lhs_str_view.data(), lhs_str_view.length(), rhs_str_view.data(),
      rhs_str_view.length(),
      edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));
}

[[nodiscard]] static auto CompressedAlignment(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads,
    biosoup::Overlap const& ovlp) -> AlignmentSummary {
  auto [lhs_substr, rhs_substr] = ExtractSubstrings(reads, ovlp);
  auto const edlib_res = AlignStrings(lhs_substr, rhs_substr);

  auto dst = AlignmentSummary{};

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

        dst.mats.push_back(MatchInterval{.first = start, .last = lhs_pos});
        break;
      }
      case 1: {
        auto start = lhs_pos;
        while (i < edlib_res.alignmentLength && edlib_res.alignment[i] == 1) {
          ++lhs_pos;
          ++i;
        }

        dst.dels.push_back(DeletinInterval{.first = start, .last = lhs_pos});
        break;
      }

      case 2: {
        auto val = 0ULL;
        auto len = 0U;
        while (i < edlib_res.alignmentLength && edlib_res.alignment[i] == 2) {
          val = len < kMaxInsLen
                    ? (((val << 2) |
                        biosoup::kNucleotideCoder[rhs_substr[rhs_pos]]))
                    : val;
          len += (len < kMaxInsLen);

          ++rhs_pos;
          ++i;
        }

        dst.inss.push_back(InsertionSite{
            .pos = lhs_pos,
            .ins =
                Insertion{.val = val, .len = static_cast<std::uint8_t>(len)}});
        break;
      }
      case 3: {
        dst.miss.push_back(MismatchSite{
            .pos = lhs_pos,
            .code = biosoup::kNucleotideCoder[rhs_substr[rhs_pos]]});
        ++i;
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
    std::vector<AlignmentSummary> const& alignments,
    std::uint32_t const kBackboneSignal)
    -> std::unique_ptr<biosoup::NucleicAcid> {
  auto covg = std::vector<CoverageSignals>(read->inflated_len);

  for (auto i = 0U; i < read->inflated_len; ++i) {
    covg[i].signals[read->Code(i)] += kBackboneSignal;
  }

  auto insertions =
      tsl::robin_map<std::uint32_t,
                     std::vector<std::pair<Insertion, std::uint32_t>>>();

  auto const record_insertion = [&insertions](InsertionSite const& is) -> void {
    auto& ins_vec = insertions[is.pos];
    auto ins_iter = std::find_if(
        ins_vec.begin(), ins_vec.end(),
        [ins = is.ins](
            std::pair<Insertion, std::uint32_t> const& ins_cnt) -> bool {
          return ins_cnt.first.val == ins.val && ins_cnt.first.len == ins.len;
        });

    if (ins_iter != ins_vec.end()) {
      ++(ins_iter->second);
    } else {
      ins_vec.emplace_back(is.ins, 1U);
    }
  };

  auto const decode_insertion =
      [&insertions](std::uint32_t const pos) -> std::string {
    auto [val, len] = insertions[pos].front().first;
    auto dst = std::string(len, '\0');

    for (auto i = 0U; i < len; ++i) {
      dst[len - (i + 1U)] = biosoup::kNucleotideDecoder[val & 3U];
      val >>= 2U;
    }

    return dst;
  };

  for (auto const& summary : alignments) {
    for (auto const [lo, hi] : summary.mats) {
      for (auto curr = lo; curr != hi; ++curr) {
        ++covg[curr].signals[read->Code(curr)];
      }
    }

    for (auto const [lo, hi] : summary.dels) {
      for (auto curr = lo; curr != hi; ++curr) {
        ++covg[curr].signals[CoverageSignals::kDelIdx];
      }
    }

    for (auto const& ins_site : summary.inss) {
      ++covg[ins_site.pos].signals[CoverageSignals::kInsIdx];
      record_insertion(ins_site);
    }

    for (auto const [pos, code] : summary.miss) {
      ++covg[pos].signals[code];
    }
  }

  // filter weak insertions
  for (auto iter = insertions.begin(); iter != insertions.end(); ++iter) {
    auto& ins_vec = iter.value();
    if (auto mx_iter = std::max_element(
            ins_vec.begin(), ins_vec.end(),
            [](std::pair<Insertion, std::uint32_t> const& lhs,
               std::pair<Insertion, std::uint32_t> const& rhs) -> bool {
              return lhs.second < rhs.second;
            });
        mx_iter != ins_vec.begin()) {
      std::swap(*mx_iter, *ins_vec.begin());
    }

    ins_vec.resize(1);
  }

  auto corrected = std::string();
  corrected.reserve(read->inflated_len + insertions.size() * 5);

  for (auto i = 0U; i < covg.size(); ++i) {
    auto const kMxIdx =
        std::distance(covg[i].signals.begin(),
                      std::max_element(covg[i].signals.begin(),
                                       covg[i].signals.begin() + 5U));

    auto const kSigSum =
        std::accumulate(covg[i].signals.cbegin(), covg[i].signals.cbegin() + 5U,
                        0U, std::plus<std::uint32_t>());

    if (kMxIdx != CoverageSignals::kDelIdx) {
      corrected.push_back(biosoup::kNucleotideDecoder[kMxIdx]);
    }

    if (covg[i].signals[CoverageSignals::kInsIdx] > kBackboneSignal &&
        covg[i].signals[CoverageSignals::kInsIdx] > 0.66 * kSigSum) {
      auto const kInsStr = decode_insertion(i);
      corrected.insert(corrected.end(), kInsStr.cbegin(), kInsStr.cend());
    }
  }

  return std::make_unique<biosoup::NucleicAcid>(read->name, corrected);
}

}  // namespace detail

auto ErrorCorrect(State& state, MapCfg const map_cfg,
                  CorrectConfig const correct_cfg,
                  std::vector<std::unique_ptr<biosoup::NucleicAcid>> src_reads)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
  auto corrected_targets = std::vector<std::unique_ptr<biosoup::NucleicAcid>>();
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
        std::vector<std::vector<detail::AlignmentSummary>>(kNTargets);

    auto align_futures = std::vector<std::future<detail::AlignmentSummary>>();

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
          [kCovgEstimate](std::unique_ptr<biosoup::NucleicAcid> const& read,
                          std::vector<detail::AlignmentSummary>& alignments)
              -> std::unique_ptr<biosoup::NucleicAcid> {
            auto dst = detail::GenerateConsensus(
                read, alignments, std::max(2U, kCovgEstimate / 3U));

            alignments.clear();
            return dst;
          },
          std::cref(src_reads[target_id]), std::ref(alignments[i])));

      if (i > 0 && (i & 127U) == 0U || i + 1U == kNTargets) {
        fmt::print(stderr,
                   "\r[camel::ErrorCorrect]({:12.3f}) collected {} / {} "
                   "alignment futures",
                   timer.Lap(), i + 1U, kNTargets);
      }
    }

    timer.Stop();
    fmt::print(stderr, "\n");

    timer.Start();
    for (auto i = 0U; i < consensus_futures.size(); ++i) {
      corrected_targets.push_back(consensus_futures[i].get());
      if (i > 0U && (i & 127U) == 0U || i + 1U == kNTargets) {
        fmt::print(stderr,
                   "\r[camel::ErrorCorrect]({:12.3f}) collected {} / {} "
                   "consensus_futures futures",
                   timer.Lap(), i + 1U, kNTargets);
      }
    }

    timer.Stop();
    fmt::print(stderr, "\n");
  }

  {
    timer.Start();
    for (auto i = 0U; i < corrected_targets.size(); ++i) {
      corrected_targets[i]->id = i;
    }

    auto queries = std::vector<std::unique_ptr<biosoup::NucleicAcid>>();

    queries.reserve(kNContained);
    std::copy_if(
        std::make_move_iterator(src_reads.begin()),
        std::make_move_iterator(src_reads.end()), std::back_inserter(queries),
        [&ovlps](std::unique_ptr<biosoup::NucleicAcid> const& acid) -> bool {
          return ovlps[acid->id].empty();
        });

    for (auto i = 0U; i < queries.size(); ++i) {
      queries[i]->id = kNTargets + i;
    }

    fmt::print(
        stderr,
        "[camel::ErrorCorrect]({:12.3f}) prepared reads for reconstruction\n",
        timer.Stop());

    timer.Start();
    auto minimizer_engine = ram::MinimizerEngine(
        state.thread_pool, map_cfg.kmer_len, map_cfg.win_len);

    minimizer_engine.Minimize(corrected_targets.cbegin(),
                              corrected_targets.cend());
    minimizer_engine.Filter(map_cfg.filter_p);

    fmt::print(
        stderr,
        "[camel::ErrorCorrect]({:12.3f}) minimized {} corrected targets\n",
        timer.Stop(), corrected_targets.size());

    timer.Start();
    auto reconstruct_futures = std::vector<std::future<void>>();
    for (auto& query : queries) {
      reconstruct_futures.emplace_back(state.thread_pool->Submit(
          [&corrected_targets, &minimizer_engine](
              std::unique_ptr<biosoup::NucleicAcid>& acid) -> void {
            auto ovlps = minimizer_engine.Map(acid, true, true, true);
            ovlps.erase(
                std::remove_if(
                    ovlps.begin(), ovlps.end(),
                    [&corrected_targets,
                     &acid](biosoup::Overlap const& ovlp) -> bool {
                      auto const kOvlpType = detail::DetermineOverlapType(
                          ovlp, corrected_targets[ovlp.rhs_id]->inflated_len,
                          acid->inflated_len);
                      return detail::OverlapError(ovlp) > 0.3 ||
                             kOvlpType != detail::OverlapType::kLhsContained;
                    }),
                ovlps.end());

            if (!ovlps.empty()) {
              auto const mx_ovlp = std::max_element(
                  ovlps.begin(), ovlps.end(),
                  [](biosoup::Overlap const& a,
                     biosoup::Overlap const& b) -> bool {
                    return detail::OverlapLength(a) < detail::OverlapLength(b);
                  });

              auto target_extract =
                  corrected_targets[mx_ovlp->rhs_id]->InflateData(
                      mx_ovlp->rhs_begin,
                      mx_ovlp->rhs_end - mx_ovlp->rhs_begin);

              if (!(mx_ovlp->strand)) {
                auto rc_buff = biosoup::NucleicAcid("", target_extract);
                rc_buff.ReverseAndComplement();

                target_extract = rc_buff.InflateData();
              }

              auto acid_corrected = std::make_unique<biosoup::NucleicAcid>(
                  acid->name, target_extract);

              acid = std::move(acid_corrected);
            } else {
              acid->block_quality.clear();  // TODO: cheeky
            }
          },
          std::ref(query)));
    }

    for (auto i = 0U; i < reconstruct_futures.size(); ++i) {
      reconstruct_futures[i].wait();
      if ((i & 1023) == 0U || i + 1U == reconstruct_futures.size()) {
        fmt::print(
            stderr,
            "\r[camel::ErrorCorrect]({:12.3f}) reconstructed {} / {} reads",
            timer.Lap(), i + 1U, reconstruct_futures.size());
      }
    }

    timer.Stop();
    fmt::print(stderr, "\n");

    src_reads.clear();

    dst.reserve(kNTargets + kNContained);
    std::move(corrected_targets.begin(), corrected_targets.end(),
              std::back_inserter(dst));
    std::move(queries.begin(), queries.end(), std::back_inserter(dst));
  }

  fmt::print(stderr, "[camel::ErrorCorrect]({:12.3f})\n", timer.elapsed_time());
  return dst;
}

}  // namespace camel
