#include "camel/correct.h"

#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <numeric>
#include <random>

#include "biosoup/timer.hpp"
#include "camel/coverage.h"
#include "camel/io.h"
#include "camel/mapping.h"
#include "detail/overlap.h"
#include "edlib.h"
#include "fmt/core.h"
#include "spoa/alignment_engine.hpp"
#include "spoa/spoa.hpp"
#include "tsl/robin_map.h"

namespace camel {

namespace detail {

static auto constexpr kWinPadding = 13U;
static auto constexpr kSpikeMergeLen = 240U;

static auto constexpr kAllowedFuzzPercent = 0.01;
static auto constexpr kSmallWindowPercent = 0.04;

static auto constexpr kReportMaks = (1U << 10U) - 1U;
static auto constexpr kSnpProximityLimit = 50U;

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
  std::uint8_t code;
};

auto operator==(Evidence const lhs, Evidence const rhs) noexcept -> bool {
  return lhs.query_pos == rhs.query_pos && lhs.code == rhs.code;
}

auto operator!=(Evidence const lhs, Evidence const rhs) noexcept -> bool {
  return !(lhs == rhs);
}

struct Segment {
  Interval query_interval;
  Interval target_interval;
  std::uint32_t target_id;

  std::uint32_t diff_score;
  std::vector<Evidence> snp_evidence;
  std::vector<Evidence> indel_evidence;
};

struct PloidyInterval : Interval {
  std::vector<Evidence> snp_evidence;
  std::vector<std::uint32_t> indel_signals;
};

struct CorrectionInterval : PloidyInterval {
  std::vector<Segment> segments;
};

struct OverlapEdlibAlignment {
  biosoup::Overlap const ovlp;
  EdlibAlignResult edlib_result;
};

struct AnnotedRead {
  std::unique_ptr<biosoup::NucleicAcid> read;
  std::vector<std::uint32_t> snp_sites;
};

static auto EvidenceDiff(std::vector<Evidence> const& lhs,
                         std::vector<Evidence> const& rhs) -> std::uint32_t {
  auto dst = 0U;
  auto i = 0U, j = 0U;
  for (; i < lhs.size() && j < rhs.size(); ++i, ++j) {
    if (lhs[i].query_pos == rhs[j].query_pos) {
      dst += lhs[i].code != rhs[j].code;
    } else if (lhs[i].query_pos < rhs[j].query_pos) {
      ++dst, ++i;
    } else {
      ++dst, ++j;
    }
  }

  dst += std::max(lhs.size() - i, rhs.size() - j);
  return dst;
}

[[nodiscard]] static inline auto IntervalLen(Interval const& intv)
    -> std::uint32_t {
  return intv.end_idx - intv.start_idx;
}

[[nodiscard]] static auto FindOverlapsAndFilterReads(
    State& state, MapCfg const map_cfg,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> src_reads)
    -> tsl::robin_map<std::uint32_t, ReadOverlapsPair> {
  auto reads_overlaps = std::vector<ReadOverlapsPair>();
  reads_overlaps.reserve(src_reads.size());

  auto overlaps = camel::FindConfidentOverlaps(state, map_cfg, src_reads);
  auto timer = biosoup::Timer();

  timer.Start();
  reads_overlaps.reserve(src_reads.size());
  std::transform(
      std::make_move_iterator(src_reads.begin()),
      std::make_move_iterator(src_reads.end()),
      std::make_move_iterator(overlaps.begin()),
      std::back_inserter(reads_overlaps),
      [](std::unique_ptr<biosoup::NucleicAcid> read,
         std::vector<biosoup::Overlap> ovlps) -> camel::ReadOverlapsPair {
        return {.read = std::move(read), .overlaps = std::move(ovlps)};
      });

  fmt::print(stderr,
             "[camel::detail::CollectOverlapsAndFilterReads]({:12.3f}) tied "
             "reads with overlaps\n",
             timer.Stop());

  {
    using namespace std::placeholders;
    using namespace std::literals;

    auto ovlp_cnsts = std::vector<std::uint64_t>(reads_overlaps.size());
    for (auto const& ro : reads_overlaps) {
      for (auto const& ovlp : ro.overlaps) {
        ++ovlp_cnsts[ovlp.lhs_id];
        ++ovlp_cnsts[ovlp.rhs_id];
      }
    }

    // preserving iterators
    auto unmapped_first = std::stable_partition(
        reads_overlaps.begin(), reads_overlaps.end(),
        [&ovlp_cnsts](camel::ReadOverlapsPair const& ro) -> bool {
          return ovlp_cnsts[ro.read->id] > 0UL;
        });

    auto const n_mapped = std::distance(reads_overlaps.begin(), unmapped_first);
    auto const n_unmapped = reads_overlaps.size() - n_mapped;

    timer.Start();

    // store unmapped reads
    if (unmapped_first != reads_overlaps.end()) {
      auto unmapped = std::vector<std::unique_ptr<biosoup::NucleicAcid>>();
      unmapped.reserve(n_unmapped);

      std::transform(std::make_move_iterator(unmapped_first),
                     std::make_move_iterator(reads_overlaps.end()),
                     std::back_inserter(unmapped),
                     [](camel::ReadOverlapsPair ro)
                         -> std::unique_ptr<biosoup::NucleicAcid> {
                       return std::move(ro.read);
                     });

      auto const dst_folder = state.log_path / "unmapped";

      if (std::filesystem::exists(dst_folder)) {
        std::filesystem::remove_all(dst_folder);
      }

      std::filesystem::create_directory(dst_folder);
      camel::StoreSequences(state, unmapped, dst_folder);

      reads_overlaps.resize(n_mapped);
    }

    fmt::print(
        stderr,
        "[camel::detail::CollectOverlapsAndFilterReads]({:12.3f}) stored {} / "
        "{} unmapped reads\n",
        timer.Stop(), n_unmapped, n_unmapped + n_mapped);
  }

  auto dst = tsl::robin_map<std::uint32_t, ReadOverlapsPair>();
  dst.reserve(reads_overlaps.size());

  for (auto& it : reads_overlaps) {
    auto id = it.read->id;
    dst[id] = std::move(it);
  }

  return dst;
}

// TODO: extract to overlap detail
[[nodiscard]] static auto OverlapStrings(
    tsl::robin_map<std::uint32_t, ReadOverlapsPair> const& reads_overlaps,
    biosoup::Overlap const& ovlp) -> std::pair<std::string, std::string> {
  auto query_str =
      reads_overlaps.at(ovlp.lhs_id)
          .read->InflateData(ovlp.lhs_begin, ovlp.lhs_end - ovlp.lhs_begin);

  auto target_str =
      reads_overlaps.at(ovlp.rhs_id)
          .read->InflateData(ovlp.rhs_begin, ovlp.rhs_end - ovlp.rhs_begin);

  if (!ovlp.strand) {
    auto rc = biosoup::NucleicAcid("", target_str);
    rc.ReverseAndComplement();

    target_str = rc.InflateData();
  }

  return std::pair(std::move(query_str), std::move(target_str));
}

[[nodiscard]] static auto AlignStrings(std::string const& query_str,
                                       std::string const& target_str)
    -> EdlibAlignResult {
  return edlibAlign(
      query_str.c_str(), query_str.size(), target_str.c_str(),
      target_str.size(),
      edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));
}

[[nodiscard]] static auto FindAlignments(
    State& state,
    tsl::robin_map<std::uint32_t, ReadOverlapsPair> const& reads_overlaps,
    std::uint32_t const query_id) -> std::vector<OverlapEdlibAlignment> {
  auto dst = std::vector<OverlapEdlibAlignment>();

  auto const create_alignment =
      [&reads_overlaps](biosoup::Overlap const& ovlp) -> OverlapEdlibAlignment {
    auto const [query_str, target_str] = OverlapStrings(reads_overlaps, ovlp);

    return {.ovlp = ovlp, .edlib_result = AlignStrings(query_str, target_str)};
  };

  dst.reserve(reads_overlaps.at(query_id).overlaps.size());
  for (auto const& ovlp : reads_overlaps.at(query_id).overlaps) {
    auto const rev_ovlp = ReverseOverlap(ovlp);
    dst.push_back(create_alignment(rev_ovlp));
  }

  for (auto const& [read_id, read_overlaps] : reads_overlaps) {
    if (read_id != query_id) {
      auto const& ovlps = read_overlaps.overlaps;
      auto first = std::lower_bound(
          ovlps.cbegin(), ovlps.cend(), query_id,

          [](biosoup::Overlap const& ovlp, std::uint32_t const lhs_id) -> bool {
            return ovlp.lhs_id < lhs_id;
          });

      for (; first != ovlps.cend() && first->lhs_id == query_id; ++first) {
        dst.push_back(create_alignment(*first));
      }
    }
  }

  return dst;
}

[[nodiscard]] static auto EstimateCoverage(
    State& state,
    tsl::robin_map<std::uint32_t, ReadOverlapsPair> const& reads_overlaps)
    -> std::uint16_t {
  auto const kSubsampleSize =
      std::max(static_cast<std::size_t>(reads_overlaps.size() * 0.10),
               std::min(100UL, reads_overlaps.size()));

  auto pile_futures = std::vector<std::future<std::uint16_t>>();
  pile_futures.reserve(kSubsampleSize);

  // generate random ids
  auto rng_ids = std::vector<std::uint32_t>(reads_overlaps.size());
  std::transform(reads_overlaps.cbegin(), reads_overlaps.cend(),
                 rng_ids.begin(),
                 [](std::pair<std::uint32_t, ReadOverlapsPair> const& ro)
                     -> std::uint32_t { return ro.first; });

  auto rng_engine = std::mt19937(42U);
  std::shuffle(rng_ids.begin(), rng_ids.end(), rng_engine);

  rng_ids.resize(kSubsampleSize);
  for (auto const rng_id : rng_ids) {
    pile_futures.emplace_back(state.thread_pool->Submit(
        [&reads_overlaps](std::uint32_t const read_id) -> std::uint16_t {
          auto constexpr kShrinkShift = 3U;

          auto const& backbone_read = reads_overlaps.at(read_id).read;
          auto pile = std::vector<std::uint16_t>(
              ((backbone_read->inflated_len) >> kShrinkShift) + 2U, 0U);

          auto events = std::vector<std::uint32_t>();
          for (auto const& ovlp : reads_overlaps.at(read_id).overlaps) {
            events.push_back(((ovlp.rhs_begin >> kShrinkShift) + 1U) << 1U);
            events.push_back((((ovlp.rhs_end >> kShrinkShift) - 1U) << 1U) |
                             1U);
          }

          for (auto const& [that_id, that_ro] : reads_overlaps) {
            if (that_id != read_id) {
              auto ovlp_iter =
                  std::lower_bound(that_ro.overlaps.cbegin(),
                                   that_ro.overlaps.cend(), backbone_read->id,
                                   [](biosoup::Overlap const& ovlp,
                                      std::uint32_t const read_id) -> bool {
                                     return ovlp.lhs_id < read_id;
                                   });

              for (; ovlp_iter != that_ro.overlaps.cend() &&
                     ovlp_iter->lhs_id == backbone_read->id;
                   ++ovlp_iter) {
                events.push_back(((ovlp_iter->lhs_begin >> kShrinkShift) + 1U)
                                 << 1U);
                events.push_back(
                    (((ovlp_iter->lhs_end >> kShrinkShift) - 1U) << 1U) | 1U);
              }
            }
          }

          auto convg = 0U;
          auto last_event = 0U;
          std::sort(events.begin(), events.end());
          for (auto const event : events) {
            for (auto i = last_event; i < (event >> 1U); ++i) {
              if (convg > 0U) {
                pile[i] =
                    std::clamp(pile[i] + convg, 0U,
                               1U * std::numeric_limits<std::uint16_t>::max());
              }
            }

            convg += 1 - (2 * (event & 1));
            last_event = (event >> 1U);
          }

          std::nth_element(pile.begin(),
                           std::next(pile.begin(), pile.size() / 2U),
                           pile.end());

          return pile[pile.size() / 2U];
        },
        rng_id));
  }

  auto medians = std::vector<std::uint16_t>(pile_futures.size());
  std::transform(pile_futures.begin(), pile_futures.end(), medians.begin(),
                 std::mem_fn(&std::future<std::uint16_t>::get));

  std::nth_element(medians.begin(),
                   std::next(medians.begin(), medians.size() / 2U),
                   medians.end());

  return medians[medians.size() / 2U];
}

[[nodiscard]] static auto FindIntervals(
    std::unique_ptr<biosoup::NucleicAcid> const& query_read,
    std::vector<OverlapEdlibAlignment> const& overlap_alignments,
    std::uint16_t kEstimatedCovg, double const kMatchRatio)
    -> std::vector<CorrectionInterval> {
  if (query_read->inflated_len < 1.2 * kSpikeMergeLen) {
    return {};
  }

  auto dst = std::vector<CorrectionInterval>();
  auto coverage = std::vector<FastCovg>(query_read->inflated_len);

  auto const kSnpCutOff = static_cast<std::uint16_t>(
      kEstimatedCovg + std::round(5 * std::sqrt(kEstimatedCovg)));

  auto const update_coverage =
      [&overlap_alignments, &coverage,
       &query_read](OverlapEdlibAlignment const& ovlp_alignment) -> void {
    auto query_pos = ovlp_alignment.ovlp.lhs_begin;
    auto const& edlib_res = ovlp_alignment.edlib_result;
    for (auto i = 0U; i < edlib_res.alignmentLength; ++i) {
      ++coverage[query_pos].signal[edlib_res.alignment[i]];
      query_pos += (edlib_res.alignment[i] != 2);
    }
  };

  std::for_each(overlap_alignments.cbegin(), overlap_alignments.cend(),
                update_coverage);

  {
    auto spikes = std::vector<std::uint32_t>();
    spikes.reserve(query_read->inflated_len / 5U);

    auto snp_candidates = std::vector<std::pair<std::uint32_t, bool>>();
    snp_candidates.reserve(query_read->inflated_len / 10U);

    auto indel_candidates = std::vector<std::uint32_t>();
    indel_candidates.reserve(snp_candidates.size());

    for (auto pos = 0U; pos < coverage.size(); ++pos) {
      auto const& covg = coverage[pos];
      if (covg.signal[0] < std::round(kMatchRatio * kEstimatedCovg)) {
        spikes.push_back(pos);
      }

      if (kEstimatedCovg * 0.25 <= covg.signal[3] &&
          covg.signal[3] <= 0.8 * kEstimatedCovg &&
          covg.signal[3] < kSnpCutOff) {
        snp_candidates.emplace_back(pos, true);
      }

      if ((kEstimatedCovg * 0.33 <= covg.signal[1] &&
           covg.signal[1] <= 0.66 * kEstimatedCovg &&
           covg.signal[1] <= kSnpCutOff) ||

          (kEstimatedCovg * 0.33 <= covg.signal[2] &&
           covg.signal[1] <= 0.66 * kEstimatedCovg &&
           covg.signal[1] <= kSnpCutOff)) {
        indel_candidates.push_back(pos);
      }
    }

    if (!spikes.empty()) {
      // create initial groups

      auto groups = std::vector<std::pair<std::uint32_t, std::uint32_t>>();
      groups.resize(spikes.size());

      std::generate_n(
          groups.begin(), groups.size(),
          [i = 0U]() mutable -> std::pair<std::uint32_t, std::uint32_t> {
            return std::pair(i, i++);
          });

      while (true) {
        auto term = true;
        for (auto i = 1U; i + 1U < groups.size(); ++i) {
          auto const lhs_dist =
              spikes[groups[i].second] - spikes[groups[i - 1U].first];
          auto const rhs_dist =
              spikes[groups[i + 1U].second] - spikes[groups[i].first];

          if (lhs_dist < rhs_dist && lhs_dist < kSpikeMergeLen) {
            groups[i].first = groups[i - 1U].first;
            groups[i - 1U] = {1U, 0U};

            term = false;
          }

          if (lhs_dist >= rhs_dist && rhs_dist < kSpikeMergeLen) {
            groups[i].second = groups[i + 1U].second;
            groups[i + 1U] = {1U, 0U};

            term = false;
            std::swap(groups[i], groups[i + 1U]);
            ++i;
          }
        }

        groups.erase(
            std::remove_if(groups.begin(), groups.end(),
                           [](std::pair<std::uint32_t, std::uint32_t> const pii)
                               -> bool { return pii.first > pii.second; }),
            groups.end());

        if (term) {
          break;
        }
      }

      dst.resize(groups.size());

      // convert groups to intervals
      dst[0].start_idx = spikes[groups[0].first] >= kWinPadding
                             ? (spikes[groups[0].first] - kWinPadding)
                             : 0U;

      for (auto i = 0U; i + 1U < groups.size(); ++i) {
        auto const rhs_gap =
            spikes[groups[i + 1U].first] - spikes[groups[i].second];

        if (rhs_gap >= 2 * kWinPadding + 1U) {
          dst[i].end_idx = spikes[groups[i].second] + kWinPadding + 1U;
          dst[i + 1U].start_idx = spikes[groups[i + 1U].first] - kWinPadding;
        } else {
          dst[i].end_idx = spikes[groups[i].second] + (rhs_gap / 2U);
          dst[i + 1U].start_idx = spikes[groups[i + 1U].first] - (rhs_gap / 2U);
        }
      }

      dst.back().end_idx = (spikes[groups.back().second] + kWinPadding <=
                            query_read->inflated_len)
                               ? (spikes[groups.back().second] + kWinPadding)
                               : query_read->inflated_len;
    }

    decltype(dst)(dst.begin(), dst.end()).swap(dst);

    // prune snps
    {
      for (auto i = 0U; i < snp_candidates.size(); ++i) {
        if (i > 0U && snp_candidates[i].second - snp_candidates[i - 1U].first <
                          kSnpProximityLimit) {
          snp_candidates[i].second = false;
        }
      }

      auto const erase_first = std::stable_partition(
          snp_candidates.begin(), snp_candidates.end(),
          [](std::pair<std::uint32_t, bool> const pub) -> bool {
            return pub.second;
          });

      snp_candidates.erase(erase_first, snp_candidates.end());
    }

    {
      auto i = 0U;
      auto j = 0U;
      for (auto& interval : dst) {
        for (; i < snp_candidates.size() &&
               snp_candidates[i].first < interval.end_idx;
             ++i) {
          if (snp_candidates[i].first >= interval.start_idx) {
            interval.snp_evidence.push_back(
                Evidence{.query_pos = snp_candidates[i].first,
                         .code = static_cast<std::uint8_t>(
                             query_read->Code(snp_candidates[i].first))});
          }
        }

        for (; j < indel_candidates.size() &&
               indel_candidates[j] < interval.end_idx;
             ++j) {
          if (indel_candidates[j] >= interval.start_idx) {
            interval.indel_signals.push_back(indel_candidates[j]);
          }
        }

        decltype(interval.snp_evidence)(interval.snp_evidence)
            .swap(interval.snp_evidence);

        decltype(interval.indel_signals)(interval.indel_signals)
            .swap(interval.indel_signals);
      }
    }
  }

  return dst;
}

[[nodiscard]] static auto CorrectRead(
    State& state,
    std::unique_ptr<spoa::AlignmentEngine> const& alignment_engine,
    tsl::robin_map<std::uint32_t, ReadOverlapsPair> const& reads_overlaps,
    std::uint32_t const query_id, std::uint16_t const kEstimatedCovg,
    double const kMatchRatio) -> std::unique_ptr<biosoup::NucleicAcid> {
  auto dst = std::unique_ptr<biosoup::NucleicAcid>();
  auto const ovlps_aligned = FindAlignments(state, reads_overlaps, query_id);

  auto const& query_read = reads_overlaps.at(query_id).read;
  auto intervals =
      FindIntervals(query_read, ovlps_aligned, kEstimatedCovg, kMatchRatio);

  if (intervals.empty()) {
    // maybe think of a better way to do this
    dst = std::make_unique<biosoup::NucleicAcid>(
        reads_overlaps.at(query_id).read->name,
        reads_overlaps.at(query_id).read->InflateData());
  } else {
    auto const backbone = query_read->InflateData();
    auto graphs = std::vector<spoa::Graph>();

    graphs.reserve(intervals.size());
    std::transform(
        intervals.cbegin(), intervals.cend(), std::back_inserter(graphs),
        [&backbone](PloidyInterval const& pi) -> spoa::Graph {
          auto graph = spoa::Graph();
          graph.AddAlignment(
              spoa::Alignment(),
              backbone.substr(pi.start_idx, pi.end_idx - pi.start_idx), 0U);

          return graph;
        });

    auto const align_to_graph = [&intervals, &graphs, &alignment_engine](
                                    std::uint32_t const interval_idx,
                                    Interval local_interval,
                                    std::string const& target_substr) -> void {
      auto const kIntervalLen = IntervalLen(intervals[interval_idx]);
      if (kIntervalLen > kWinPadding &&
          IntervalLen(local_interval) < kIntervalLen * kSmallWindowPercent) {
        return;
      }

      auto const kLegalStart = kIntervalLen * kAllowedFuzzPercent;
      auto const kLegalEnd = kIntervalLen - kLegalStart;

      auto alignment = spoa::Alignment();
      if (local_interval.start_idx <= kLegalStart && kLegalEnd <= kLegalEnd) {
        alignment =
            alignment_engine->Align(target_substr, graphs[interval_idx]);
      } else {
        auto mapping = std::vector<spoa::Graph::Node const*>();
        auto subgraph = graphs[interval_idx].Subgraph(
            local_interval.start_idx, local_interval.end_idx - 1, &mapping);
        alignment = alignment_engine->Align(target_substr, subgraph);
        subgraph.UpdateAlignment(mapping, &alignment);
      }

      graphs[interval_idx].AddAlignment(alignment, target_substr);
    };

    for (auto const& [ovlp, edlib_res] : ovlps_aligned) {
      auto intv_iter = std::upper_bound(
          intervals.begin(), intervals.end(), ovlp.lhs_begin,
          [](std::uint32_t const pos, Interval const& intv) -> bool {
            return pos < intv.end_idx;
          });

      auto snp_idx = 0U;
      auto snp_evidence = std::vector<Evidence>();

      auto indle_idx = 0U;
      auto indle_evidence = std::vector<Evidence>();

      if (intv_iter != intervals.end()) {
        snp_evidence.reserve(intv_iter->snp_evidence.size());
        indle_evidence.reserve(intv_iter->indel_signals.size());
      }

      auto query_anchor = ovlp.lhs_begin;
      auto target_anchor = ovlp.rhs_begin;

      auto query_pos = ovlp.lhs_begin;
      auto target_pos = ovlp.rhs_begin;

      auto const& target_read = reads_overlaps.at(ovlp.rhs_id).read;

      for (auto i = 0U;
           i < edlib_res.alignmentLength && intv_iter != intervals.end(); ++i) {
        if (snp_idx < intv_iter->snp_evidence.size() &&
            intv_iter->snp_evidence[snp_idx].query_pos == query_pos) {
          if (ovlp.strand) {
            snp_evidence.push_back(Evidence{
                .query_pos = query_pos,
                .code =
                    static_cast<std::uint8_t>(target_read->Code(target_pos))});
          } else {
            snp_evidence.push_back(
                Evidence{.query_pos = query_pos,
                         .code = static_cast<std::uint8_t>(
                             3 ^ target_read->Code(target_read->inflated_len -
                                                   1U - target_pos))});
          }

          ++snp_idx;
        }

        if (indle_idx < intv_iter->indel_signals.size() &&
            intv_iter->indel_signals[indle_idx] == query_pos) {
          if (ovlp.strand) {
            indle_evidence.push_back(Evidence{
                .query_pos = query_pos,
                .code =
                    static_cast<std::uint8_t>(target_read->Code(target_pos))});
          } else {
            indle_evidence.push_back(
                Evidence{.query_pos = query_pos,
                         .code = static_cast<std::uint8_t>(
                             3 ^ target_read->Code(target_read->inflated_len -
                                                   1U - target_pos))});
          }

          ++indle_idx;
        }

        query_pos += (edlib_res.alignment[i] != 2U);
        target_pos += (edlib_res.alignment[i] != 1U);

        if (query_pos == intv_iter->end_idx) {
          intv_iter->segments.push_back(
              {.query_interval = {.start_idx =
                                      query_anchor - intv_iter->start_idx,
                                  .end_idx = query_pos - intv_iter->start_idx},
               .target_interval = {.start_idx = target_anchor,
                                   .end_idx = target_pos},
               .target_id = ovlp.rhs_id,
               .snp_evidence = std::move(snp_evidence),
               .indel_evidence = std::move(indle_evidence)});

          ++intv_iter;
        }

        if (intv_iter != intervals.end() && query_pos == intv_iter->start_idx) {
          snp_idx = 0U;
          snp_evidence.reserve(intv_iter->snp_evidence.size());

          indle_idx = 0U;
          indle_evidence.reserve(intv_iter->indel_signals.size());

          query_anchor = query_pos;
          target_anchor = target_pos;
        }
      }
    }

    for (auto idx = 0U; idx < intervals.size(); ++idx) {
      for (auto& segment : intervals[idx].segments) {
        segment.diff_score =
            EvidenceDiff(segment.snp_evidence, intervals[idx].snp_evidence);
      }

      {
        auto const pivot = std::next(intervals[idx].segments.begin(),
                                     intervals[idx].segments.size() * 0.5);
        if (pivot != intervals[idx].segments.begin() &&
            pivot != intervals[idx].segments.end()) {
          std::partial_sort(intervals[idx].segments.begin(), pivot,
                            intervals[idx].segments.end(),
                            [](Segment const& lhs, Segment const& rhs) -> bool {
                              return lhs.diff_score < rhs.diff_score;
                            });

          intervals[idx].segments.erase(pivot, intervals[idx].segments.end());
        }
      }

      if (intervals[idx].segments.empty()) {
        auto const& ref_indles = intervals[idx].segments.front().indel_evidence;
        for (auto& segment : intervals[idx].segments) {
          segment.diff_score = EvidenceDiff(segment.indel_evidence, ref_indles);
        }

        auto const pivot = std::next(intervals[idx].segments.begin(),
                                     intervals[idx].segments.size() * 0.7);
        if (pivot != intervals[idx].segments.begin() &&
            pivot != intervals[idx].segments.end()) {
          std::partial_sort(intervals[idx].segments.begin(), pivot,
                            intervals[idx].segments.end(),
                            [](Segment const& lhs, Segment const& rhs) -> bool {
                              return lhs.diff_score < rhs.diff_score;
                            });

          intervals[idx].segments.erase(pivot, intervals[idx].segments.end());
        }
      }

      for (auto i = 0U; i < intervals[idx].segments.size(); ++i) {
        auto const local_interval = intervals[idx].segments[i].query_interval;
        auto const target_interval = intervals[idx].segments[i].target_interval;
        auto const target_substr =
            reads_overlaps.at(intervals[idx].segments[i].target_id)
                .read->InflateData(target_interval.start_idx,
                                   IntervalLen(target_interval));

        align_to_graph(idx, local_interval, target_substr);
      }
    }

    // generate consensus
    {
      auto consensus = std::string();
      consensus.reserve(query_read->inflated_len * 1.1);

      auto pos = 0U;
      for (auto i = 0U; i < intervals.size(); ++i) {
        auto const& intv = intervals[i];
        consensus += query_read->InflateData(pos, intv.start_idx - pos);
        consensus += graphs[i].GenerateConsensus();

        pos = intv.end_idx;
      }

      consensus += query_read->InflateData(pos);

      dst = std::make_unique<biosoup::NucleicAcid>(query_read->name, consensus);
    }
  }

  for (auto const& [ovlp, edlib_res] : ovlps_aligned) {
    edlibFreeAlignResult(edlib_res);
  }

  return dst;
}

}  // namespace detail

auto SnpErrorCorrect(
    State& state, MapCfg const map_cfg, PolishConfig const polish_cfg,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> src_reads)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
  auto dst = std::vector<std::unique_ptr<biosoup::NucleicAcid>>();

  auto timer = biosoup::Timer();

  timer.Start();
  auto reads_overlaps =
      detail::FindOverlapsAndFilterReads(state, map_cfg, std::move(src_reads));
  fmt::print(stderr,
             "[camel::SnpErrorCorrect]({:12.3f}) tied reads with overlaps\n",
             timer.Stop());

  timer.Start();
  auto const kEstimatedCovg = detail::EstimateCoverage(state, reads_overlaps);
  fmt::print(stderr,
             "[camel::SnpErrorCorrect]({:12.3f}) estimated coverage: {}\n",
             timer.Stop(), kEstimatedCovg);

  {
    timer.Start();

    auto alignment_engines =
        tsl::robin_map<std::thread::id,
                       std::unique_ptr<spoa::AlignmentEngine>>();

    for (auto const& [thread_id, _] : state.thread_pool->thread_map()) {
      alignment_engines[thread_id] = spoa::AlignmentEngine::Create(
          spoa::AlignmentType::kNW, polish_cfg.poa_cfg.match,
          polish_cfg.poa_cfg.mismatch, polish_cfg.poa_cfg.gap);
    }

    auto const find_batch_last =
        [](std::vector<ReadOverlapsPair>::const_iterator first,
           std::vector<ReadOverlapsPair>::const_iterator last,
           std::size_t const batch_cap)
        -> std::vector<ReadOverlapsPair>::const_iterator {
      for (auto batch_sz = 0UL; batch_sz < batch_cap && first != last;
           ++first) {
        batch_sz += first->read->inflated_len;
      }

      return first;
    };

    auto correct_futures =
        std::vector<std::future<std::unique_ptr<biosoup::NucleicAcid>>>();

    auto const correct_async = [&state, &reads_overlaps, &alignment_engines,
                                kEstimatedCovg](std::uint32_t const query_id)
        -> std::future<std::unique_ptr<biosoup::NucleicAcid>> {
      return state.thread_pool->Submit(
          [&alignment_engines, kEstimatedCovg](
              State& state,
              tsl::robin_map<std::uint32_t, ReadOverlapsPair> const&
                  reads_overlaps,
              std::uint32_t const& query_idx)
              -> std::unique_ptr<biosoup::NucleicAcid> {
            auto const& alignment_engine =
                alignment_engines[std::this_thread::get_id()];
            return detail::CorrectRead(state, alignment_engine, reads_overlaps,
                                       query_idx, kEstimatedCovg, 0.95);
          },
          std::ref(state), std::cref(reads_overlaps), query_id

      );
    };

    dst.reserve(reads_overlaps.size());
    correct_futures.reserve(reads_overlaps.size());
    for (auto const& [id, _] : reads_overlaps) {
      correct_futures.emplace_back(correct_async(id));
    }

    for (auto idx = 0U; idx < reads_overlaps.size(); ++idx) {
      dst.push_back(correct_futures[idx].get());

      if (((idx + 1U) & detail::kReportMaks) == 0U ||
          idx + 1U == reads_overlaps.size()) {
        fmt::print(
            stderr,
            "\r[camel::SnpErrorCorrect]({:12.3f}) corrected {} / {} reads",
            timer.Lap(), idx + 1U, reads_overlaps.size());
      }
    }

    fmt::print(stderr, "\n");
    timer.Stop();
  }

  fmt::print(stderr, "[camel::SnpErrorCorrect]({:12.3f})\n",
             timer.elapsed_time());

  return dst;
}

}  // namespace camel
