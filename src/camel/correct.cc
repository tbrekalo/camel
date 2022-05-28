#include "camel/correct.h"

#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <numeric>

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

static auto constexpr kWinPadding = 19U;
static auto constexpr kSpikeMergeLen = 420U;

static auto constexpr kAllowedFuzzPercent = 0.01;
static auto constexpr kSmallWindowPercent = 0.04;

static auto constexpr kBatchCap = 1UL << 28UL;

static auto constexpr kQuerySectionCap = 1U << 4U;

static auto constexpr kReportMaks = (1U << 10U) - 1U;

struct FastCovg {
  // mat, del, ins, mis
  std::array<std::uint_fast16_t, 4> signal;
};

struct Interval {
  std::uint32_t start_idx;
  std::uint32_t end_idx;
};

struct Section {
  Interval query_interval;
  Interval target_interval;
  std::uint32_t target_id;
};

struct PloidyInterval : Interval {
  std::vector<std::uint32_t> snp_sites;
};

struct CorrectionInterval : PloidyInterval {
  std::array<Section, kQuerySectionCap + 1U> sections;
  std::uint32_t n_active;
};

struct OverlapEdlibAlignment {
  biosoup::Overlap const ovlp;
  EdlibAlignResult edlib_result;
};

static inline auto IntervalLen(Interval const& intv) -> std::uint32_t {
  return intv.end_idx - intv.start_idx;
}

struct Alignment {
  std::uint32_t query_start;
  std::string target_str;
  EdlibAlignResult edlib_res;
};

[[nodiscard]] static auto FindOverlapsAndFilterReads(
    State& state, MapCfg const map_cfg,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> src_reads)
    -> std::vector<ReadOverlapsPair> {
  auto reads_overlaps = std::vector<ReadOverlapsPair>();
  reads_overlaps.reserve(src_reads.size());

  auto overlaps = camel::FindOverlaps(state, map_cfg, src_reads);
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

      decltype(reads_overlaps)(std::make_move_iterator(reads_overlaps.begin()),
                               std::make_move_iterator(unmapped_first))
          .swap(reads_overlaps);
    }

    fmt::print(stderr,
               "[camel::CollectOverlapsAndFilterReads]({:12.3f}) stored {} / "
               "{} unmapped reads\n",
               timer.Stop(), n_unmapped, n_unmapped + n_mapped);
  }

  return reads_overlaps;
}

// TODO: extract to overlap detail
[[nodiscard]] static auto OverlapStrings(
    std::vector<ReadOverlapsPair> const& reads_overlaps,
    biosoup::Overlap const& ovlp) -> std::pair<std::string, std::string> {
  auto const query_read_iter = std::lower_bound(
      reads_overlaps.cbegin(), reads_overlaps.cend(), ovlp.lhs_id,
      [](ReadOverlapsPair const& ro, std::uint32_t const query_id) -> bool {
        return ro.read->id < query_id;
      });

  auto const target_read_iter = std::lower_bound(
      reads_overlaps.cbegin(), reads_overlaps.cend(), ovlp.rhs_id,
      [](ReadOverlapsPair const& ro, std::uint32_t const target_id) -> bool {
        return ro.read->id < target_id;
      });

  auto query_str = query_read_iter->read->InflateData(
      ovlp.lhs_begin, ovlp.lhs_end - ovlp.lhs_begin);

  auto target_str = target_read_iter->read->InflateData(
      ovlp.rhs_begin, ovlp.rhs_end - ovlp.rhs_begin);

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
    State& state, std::vector<ReadOverlapsPair> const& reads_overlaps,
    std::uint32_t const query_read_idx) -> std::vector<OverlapEdlibAlignment> {
  auto dst = std::vector<OverlapEdlibAlignment>();
  auto const query_id = reads_overlaps[query_read_idx].read->id;

  auto const create_alignment =
      [&reads_overlaps](biosoup::Overlap const& ovlp) -> OverlapEdlibAlignment {
    auto const [query_str, target_str] = OverlapStrings(reads_overlaps, ovlp);

    return {.ovlp = ovlp, .edlib_result = AlignStrings(query_str, target_str)};
  };

  dst.reserve(reads_overlaps[query_read_idx].overlaps.size());
  for (auto const& ovlp : reads_overlaps[query_read_idx].overlaps) {
    auto const rev_ovlp = ReverseOverlap(ovlp);
    dst.push_back(create_alignment(rev_ovlp));
  }

  for (auto read_idx = 0U; read_idx < reads_overlaps.size(); ++read_idx) {
    if (read_idx != query_read_idx) {
      auto const& ovlps = reads_overlaps[read_idx].overlaps;
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

[[nodiscard]] static auto FindIntervals(
    std::unique_ptr<biosoup::NucleicAcid> const& query_read,
    std::vector<OverlapEdlibAlignment> const& overlap_alignments,
    double const match_ratio) -> std::vector<CorrectionInterval> {
  if (query_read->inflated_len < 2 * kSpikeMergeLen) {
    return {};
  }

  auto dst = std::vector<CorrectionInterval>();
  auto coverage = std::vector<FastCovg>(query_read->inflated_len);

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

  {
    auto spikes = std::vector<std::uint32_t>();
    spikes.reserve(query_read->inflated_len / 5U);

    auto snp_candidates = std::vector<std::uint32_t>();
    snp_candidates.reserve(query_read->inflated_len / 10U);

    for (auto pos = 0U; pos < coverage.size(); ++pos) {
      auto const& covg = coverage[pos];
      auto const sum = std::accumulate(covg.signal.cbegin(), covg.signal.cend(),
                                       std::uint_fast16_t(0));
      if (covg.signal[0] <
          static_cast<std::uint32_t>(std::round(match_ratio * sum))) {
        spikes.push_back(pos);
      }

      if (std::abs(static_cast<std::int32_t>(covg.signal[0]) -
                   static_cast<std::int32_t>(covg.signal[3])) <
              static_cast<std::uint32_t>(std::round(sum * 0.333)) &&
          covg.signal[0] + covg.signal[3] >
              static_cast<std::uint32_t>(std::round(sum * 0.8))) {
        snp_candidates.push_back(pos);
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

    {
      auto i = 0U;
      for (auto& interval : dst) {
        for (;
             i < snp_candidates.size() && snp_candidates[i] < interval.end_idx;
             ++i) {
          if (snp_candidates[i] >= interval.start_idx) {
            interval.snp_sites.push_back(snp_candidates[i]);
          }
        }

        decltype(interval.snp_sites)(interval.snp_sites)
            .swap(interval.snp_sites);
      }
    }
  }

  return dst;
}

[[nodiscard]] static auto CorrectRead(
    State& state,
    std::unique_ptr<spoa::AlignmentEngine> const& alignment_engine,
    std::vector<ReadOverlapsPair> const& reads_overlaps,
    std::uint32_t const query_idx) -> std::unique_ptr<biosoup::NucleicAcid> {
  auto dst = std::unique_ptr<biosoup::NucleicAcid>();
  auto const ovlps_aligned = FindAlignments(state, reads_overlaps, query_idx);

  auto const& query_read = reads_overlaps[query_idx].read;
  auto intervals = FindIntervals(query_read, ovlps_aligned, 0.9);

  if (intervals.empty()) {
    // maybe think of a better way to do this
    dst = std::make_unique<biosoup::NucleicAcid>(
        reads_overlaps[query_idx].read->name,
        reads_overlaps[query_idx].read->InflateData());
  } else {
    auto const backbone = reads_overlaps[query_idx].read->InflateData();
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

    auto unfilled = intervals.size();
    for (auto const& [ovlp, edlib_res] : ovlps_aligned) {
      auto intv_iter = std::upper_bound(
          intervals.begin(), intervals.end(), ovlp.lhs_begin,
          [](std::uint32_t const pos, Interval const& intv) -> bool {
            return pos < intv.end_idx;
          });

      auto valid = true;
      auto snp_idx = 0U;

      auto query_pos = ovlp.lhs_begin;
      auto target_pos = ovlp.rhs_begin;

      auto query_interval = Interval{query_pos, query_pos};
      auto target_interval = Interval{target_pos, target_pos};

      for (auto i = 0U;
           i < edlib_res.alignmentLength && intv_iter != intervals.end(); ++i) {
        query_pos += (edlib_res.alignment[i] != 2U);
        target_pos += (edlib_res.alignment[i] != 1U);

        valid &= (snp_idx >= intv_iter->snp_sites.size() ||
                  !(intv_iter->snp_sites[snp_idx] == query_pos &&
                    edlib_res.alignment[i] != 0U));

        snp_idx += (snp_idx < intv_iter->snp_sites.size() &&
                    intv_iter->snp_sites[snp_idx] == query_pos);

        if (query_pos == intv_iter->end_idx) {
          if (valid) {
            intv_iter->sections[intv_iter->n_active] = {
                .query_interval = query_interval,
                .target_interval = target_interval,
                .target_id = ovlp.rhs_id};

            unfilled -= (intv_iter->n_active + 1U == kQuerySectionCap);
            intv_iter->n_active += (intv_iter->n_active != kQuerySectionCap);
          }

          ++intv_iter;
        }

        if (query_pos == intv_iter->start_idx) {
          valid = true;
          snp_idx = 0U;

          query_interval = {query_pos, query_pos};
          target_interval = {target_pos, target_pos};
        }
      }

      if (unfilled == 0U) {
        break;
      }
    }

    for (auto idx = 0U; idx < intervals.size(); ++idx) {
      for (auto i = 0U; i < intervals[idx].n_active; ++i) {
        auto const local_interval = intervals[idx].sections[idx].query_interval;
        auto const target_interval =
            intervals[idx].sections[idx].target_interval;
        auto const target_substr =
            std::lower_bound(reads_overlaps.begin(), reads_overlaps.end(),
                             intervals[idx].sections[i].target_id,
                             [](ReadOverlapsPair const& ro,
                                std::uint32_t const target_id) -> bool {
                               return ro.read->id < target_id;
                             })
                ->read->InflateData(target_interval.start_idx,
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

    auto const correct_async = [&state, &reads_overlaps, &alignment_engines](
                                   std::uint32_t const query_idx)
        -> std::future<std::unique_ptr<biosoup::NucleicAcid>> {
      return state.thread_pool->Submit(
          [&alignment_engines](
              State& state, std::vector<ReadOverlapsPair> const& reads_overlaps,
              std::uint32_t const& query_idx)
              -> std::unique_ptr<biosoup::NucleicAcid> {
            auto const& alignment_engine =
                alignment_engines[std::this_thread::get_id()];
            return detail::CorrectRead(state, alignment_engine, reads_overlaps,
                                       query_idx);
          },
          std::ref(state), std::cref(reads_overlaps), query_idx

      );
    };

    dst.reserve(reads_overlaps.size());
    correct_futures.reserve(reads_overlaps.size());
    for (auto idx = 0U; idx < reads_overlaps.size(); ++idx) {
      correct_futures.emplace_back(correct_async(idx));
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
