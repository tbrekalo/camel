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

static auto constexpr kPackSize = 100U;

static auto constexpr kWinPadding = 19U;
static auto constexpr kSpikeMergeLen = 420U;

static auto constexpr kAllowedFuzzPercent = 0.01;
static auto constexpr kSmallWindowPercent = 0.04;

static auto constexpr kBatchCap = 1UL << 24UL;

struct FastCovg {
  // mat, del, ins, mis
  std::array<std::uint_fast16_t, 4> signal;
};

struct Interval {
  std::uint32_t start_idx;
  std::uint32_t end_idx;
};

struct PloidyInterval : Interval {
  std::vector<std::uint32_t> snp_sites;
};

static auto IntervalLen(Interval const& intv) -> std::uint32_t {
  return intv.end_idx - intv.start_idx;
}

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

[[nodiscard]] static auto FindIntervals(
    State& state, std::vector<ReadOverlapsPair> const& reads_overlaps,
    std::uint32_t query_read_idx, double const match_ratio)
    -> std::vector<PloidyInterval> {
  auto const& query_read = reads_overlaps[query_read_idx].read;
  if (query_read->inflated_len < 2 * kSpikeMergeLen) {
    return {};
  }

  auto dst = std::vector<PloidyInterval>();
  auto coverage = std::vector<FastCovg>(query_read->inflated_len);

  auto const update_coverage = [&reads_overlaps, &coverage, &query_read](
                                   biosoup::Overlap const& ovlp) -> void {
    auto const [query_substr, target_substr] =
        OverlapStrings(reads_overlaps, ovlp);
    auto const edlib_res = AlignStrings(query_substr, target_substr);

    auto query_pos = ovlp.lhs_begin;
    for (auto i = 0U; i < edlib_res.alignmentLength; ++i) {
      ++coverage[query_pos].signal[edlib_res.alignment[i]];
      query_pos += (edlib_res.alignment[i] != 2);
    }

    edlibFreeAlignResult(edlib_res);
  };

  for (auto const& [target_read, target_overlaps] : reads_overlaps) {
    if (target_read->id == query_read->id) {
      for (auto const& ovlp : target_overlaps) {
        update_coverage(detail::ReverseOverlap(ovlp));
      }
    } else {
      auto ovlp_iter = std::lower_bound(
          target_overlaps.cbegin(), target_overlaps.cend(), target_read->id,
          [](biosoup::Overlap const& ovlp, std::uint32_t const read_id)
              -> bool { return ovlp.lhs_id < read_id; });

      for (; ovlp_iter != target_overlaps.cend() &&
             ovlp_iter->lhs_id == target_read->id;
           ++ovlp_iter) {
        update_coverage(*ovlp_iter);
      }
    }
  }

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
    std::uint32_t const seq_idx) -> AnnotatedRead {
  auto dst = AnnotatedRead();
  auto const intervals = FindIntervals(state, reads_overlaps, seq_idx, 0.9);

  if (intervals.empty()) {
    // maybe think of a better way to do this
    dst.read = std::make_unique<biosoup::NucleicAcid>(
        reads_overlaps[seq_idx].read->name,
        reads_overlaps[seq_idx].read->InflateData());
  }

  auto const& query_read = reads_overlaps[seq_idx].read;

  auto const backbone = reads_overlaps[seq_idx].read->InflateData();
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

  // auto const is_lhs_usable_ovlp =
  //     [&intervals](biosoup::Overlap const& ovlp) -> bool {
  //   auto interval_iter = std::lower_bound(
  //       intervals.cbegin(), intervals.cend(), ovlp.lhs_end,
  //       [](PloidyInterval const& pi, std::uint32_t const pos) -> bool {
  //         return pi.start_idx < pos;
  //       });

  //   return interval_iter != intervals.cend() &&
  //          ovlp.lhs_begin < interval_iter->end_idx;
  // };

  // auto const is_rhs_usable_ovlp =
  //     [&intervals](biosoup::Overlap const& ovlp) -> bool {
  //   auto interval_iter = std::lower_bound(
  //       intervals.cbegin(), intervals.cend(), ovlp.rhs_end,
  //       [](PloidyInterval const& pi, std::uint32_t const pos) -> bool {
  //         return pi.start_idx < pos;
  //       });

  //   return interval_iter != intervals.cend() &&
  //          ovlp.rhs_begin < interval_iter->end_idx;
  // };

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
      alignment = alignment_engine->Align(target_substr, graphs[interval_idx]);
    } else {
      auto mapping = std::vector<spoa::Graph::Node const*>();
      auto subgraph = graphs[interval_idx].Subgraph(
          local_interval.start_idx, local_interval.end_idx - 1, &mapping);
      alignment = alignment_engine->Align(target_substr, subgraph);
      subgraph.UpdateAlignment(mapping, &alignment);
    }

    graphs[interval_idx].AddAlignment(alignment, target_substr);
  };

  auto const align_to_intervals =
      [&reads_overlaps, &intervals, &graphs,
       &align_to_graph](biosoup::Overlap const& ovlp) -> void {
    auto const [query_str, target_str] = OverlapStrings(reads_overlaps, ovlp);

    auto const edlib_res = AlignStrings(query_str, target_str);

    auto query_pos = ovlp.lhs_begin;
    auto target_pos = 0U;

    auto query_anchor = query_pos;
    auto target_anchor = target_pos;

    auto interval_idx = std::distance(
        intervals.cbegin(),
        std::upper_bound(intervals.cbegin(), intervals.cend(), query_pos,
                         [](std::uint32_t const pos, PloidyInterval const& pi)
                             -> bool { return pos < pi.end_idx; }));

    for (auto i = 0U;
         i < edlib_res.alignmentLength && interval_idx < intervals.size();
         ++i) {
      query_pos += (edlib_res.alignment[i] != 2);
      target_pos += (edlib_res.alignment[i] != 1);

      if (intervals[interval_idx].end_idx == query_pos) {
        auto const local_interval = Interval{
            .start_idx = query_anchor - intervals[interval_idx].start_idx,
            .end_idx = query_pos - intervals[interval_idx].start_idx};

        auto const target_substr =
            target_str.substr(target_anchor, target_pos - target_anchor);

        align_to_graph(interval_idx, local_interval, target_substr);
        ++interval_idx;
      }

      if (interval_idx != intervals.size() &&
          intervals[interval_idx].start_idx == query_pos) {
        query_anchor = query_pos;
        target_anchor = target_pos;
      }
    }
  };

  for (auto const& ovlp : reads_overlaps[seq_idx].overlaps) {
    auto const rev_ovlp = detail::ReverseOverlap(ovlp);
    align_to_intervals(rev_ovlp);
  }

  for (auto idx = 0U; idx != reads_overlaps.size(); ++idx) {
    if (idx != seq_idx) {
      auto curr_ovlp_iter = std::lower_bound(
          reads_overlaps[idx].overlaps.cbegin(),
          reads_overlaps[idx].overlaps.cend(), query_read->id,
          [](biosoup::Overlap const& ovlp, std::uint32_t const lhs_id) -> bool {
            return ovlp.lhs_id < lhs_id;
          });

      for (; curr_ovlp_iter != reads_overlaps[idx].overlaps.cend() &&
             curr_ovlp_iter->lhs_id == query_read->id;
           ++curr_ovlp_iter) {
        align_to_intervals(*curr_ovlp_iter);
      }
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

    dst.read =
        std::make_unique<biosoup::NucleicAcid>(query_read->name, consensus);
  }

  return dst;
}

}  // namespace detail

auto SnpErrorCorrect(
    State& state, MapCfg const map_cfg, PolishConfig const polish_cfg,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> src_reads)
    -> std::vector<AnnotatedRead> {
  auto dst = std::vector<AnnotatedRead>();

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

    auto worker_futures = std::vector<std::future<void>>();
    auto atomic_idx = std::atomic<std::uint32_t>{0U};

    auto mtx = std::mutex();
    auto report_cv = std::condition_variable();

    auto report_val = 0U;

    fmt::print(stderr, "[camel::SnpErrorCorrect] is indexing lock free: {}\n",
               decltype(atomic_idx)::is_always_lock_free);

    timer.Start();
    dst.resize(reads_overlaps.size());
    worker_futures.reserve(state.thread_pool->num_threads());
    std::generate_n(
        std::back_inserter(worker_futures), state.thread_pool->num_threads(),
        [&]() -> std::future<void> {
          return state.thread_pool->Submit([&]() -> void {
            auto const& align_engine =
                alignment_engines[std::this_thread::get_id()];

            auto curr_idx = 0U;
            while (true) {
              curr_idx = atomic_idx.fetch_add(1);
              if (curr_idx >= reads_overlaps.size()) {
                break;
              }

              dst[curr_idx] = detail::CorrectRead(state, align_engine,
                                                  reads_overlaps, curr_idx);

              if (curr_idx > 0 && ((curr_idx + 1U) % 1000 == 0U ||
                                   curr_idx + 1U == reads_overlaps.size())) {
                auto lk = std::unique_lock(mtx);
                report_val = curr_idx + 1U;
                report_cv.notify_one();
              }
            }
          });
        });

    while (true) {
      auto lk = std::unique_lock(mtx);
      report_cv.wait(lk);

      fmt::print(stderr,
                 "\r[camel::SnpErrorCorrect]({:12.3f}) corrected {} / {} reads",
                 timer.Lap(), report_val, reads_overlaps.size());

      if (report_val == reads_overlaps.size()) {
        fmt::print(stderr, "\n");
        timer.Stop();
        break;
      }
    }

    std::for_each(worker_futures.begin(), worker_futures.end(),
                  std::mem_fn(&std::future<void>::wait));

    timer.Stop();
    fmt::print(stderr, "\n");
  }

  fmt::print(stderr, "[camel::SnpErrorCorrect]({:12.3f})\n",
             timer.elapsed_time());

  return dst;
}

}  // namespace camel
