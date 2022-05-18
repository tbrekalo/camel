#include "camel/correct.h"

#include <array>
#include <chrono>
#include <cmath>
#include <deque>
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

namespace camel {

namespace detail {

static auto constexpr kPackSize = 100U;

static auto constexpr kWinPadding = 19U;
static auto constexpr kSpikeMergeLen = 420;

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

[[nodiscard]] static auto FindOverlapsAndFilterReads(
    State& state, std::vector<std::unique_ptr<biosoup::NucleicAcid>> src_reads)
    -> std::vector<ReadOverlapsPair> {
  auto reads_overlaps = std::vector<ReadOverlapsPair>();
  reads_overlaps.reserve(src_reads.size());

  auto overlaps = camel::FindOverlaps(state, MapCfg{}, src_reads);
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

struct CmpOvlpLhsReadId {
  auto operator()(biosoup::Overlap const& ovlp,
                  std::uint32_t const read_id) const noexcept -> bool {
    return ovlp.lhs_id < read_id;
  }

  auto operator()(std::uint32_t const read_id,
                  biosoup::Overlap const& ovlp) const noexcept -> bool {
    return read_id < ovlp.lhs_id;
  }
};

[[nodiscard]] static auto CollectOverlaps(
    State& state, std::vector<ReadOverlapsPair> const& ros,
    std::uint32_t const read_id) -> std::vector<biosoup::Overlap> {
  static constexpr auto kIntervalLen = 5000U;

  auto dst = std::vector<biosoup::Overlap>();
  auto const kNIntervals = std::ceil((1.0 * ros.size()) / kIntervalLen);

  auto interval_counts = std::vector<std::uint64_t>(kNIntervals);

  {
    auto count_futures = std::vector<std::future<std::uint64_t>>();

    count_futures.reserve(kNIntervals);
    for (auto first = ros.cbegin(); first != ros.cend();) {
      auto const last = std::min(std::next(first, kIntervalLen), ros.cend());

      count_futures.emplace_back(state.thread_pool->Submit(
          [read_id](std::vector<ReadOverlapsPair>::const_iterator first,
                    std::vector<ReadOverlapsPair>::const_iterator last)
              -> std::uint64_t {
            return std::transform_reduce(
                first, last, 0UL, std::plus<std::uint64_t>(),
                [read_id](ReadOverlapsPair const& ro) -> std::uint64_t {
                  auto const& ovlps = ro.overlaps;
                  if (!ovlps.empty() && ovlps.front().rhs_id == read_id) {
                    return ovlps.size();
                  } else {
                    auto const [lo, hi] =
                        std::equal_range(ovlps.cbegin(), ovlps.cend(), read_id,
                                         CmpOvlpLhsReadId());

                    return std::distance(lo, hi);
                  }
                });
          },
          first, last));

      first = last;
    }

    std::transform(
        count_futures.begin(), count_futures.end(), interval_counts.begin(),
        [](std::future<std::uint64_t>& f) -> std::uint64_t { return f.get(); });

    auto const n_total_ovlps =
        std::accumulate(interval_counts.cbegin(), interval_counts.cend(), 0UL);

    dst.resize(n_total_ovlps);
  }

  {
    auto collect_futures = std::vector<std::future<void>>();
    collect_futures.reserve(kNIntervals);

    auto interval_firsts = std::move(interval_counts);
    std::exclusive_scan(interval_firsts.begin(), interval_firsts.end(),
                        interval_firsts.begin(), 0UL);

    auto interval_id = 0U;
    for (auto first = ros.cbegin(); first != ros.cend(); ++interval_id) {
      auto const last = std::min(std::next(first, kIntervalLen), ros.cend());

      collect_futures.emplace_back(state.thread_pool->Submit(
          [read_id](std::vector<ReadOverlapsPair>::const_iterator first,
                    std::vector<ReadOverlapsPair>::const_iterator last,
                    std::vector<biosoup::Overlap>::iterator first_out) mutable
          -> void {
            for (; first != last; ++first) {
              auto const& ovlps = first->overlaps;
              if (!ovlps.empty() && ovlps.front().rhs_id == read_id) {
                std::transform(ovlps.begin(), ovlps.end(), first_out,
                               detail::ReverseOverlap);
                std::advance(first_out, ovlps.size());

              } else {
                auto const lo =
                    std::lower_bound(ovlps.cbegin(), ovlps.cend(), read_id,
                                     [](biosoup::Overlap const& ovlp,
                                        std::uint32_t const read_id) -> bool {
                                       return ovlp.lhs_id < read_id;
                                     });

                auto const hi = std::find_if_not(
                    lo, ovlps.cend(),
                    [read_id](biosoup::Overlap const& ovlp) -> bool {
                      return ovlp.lhs_id == read_id;
                    });

                first_out = std::copy(lo, hi, first_out);
              }
            }
          },
          first, last, std::next(dst.begin(), interval_firsts[interval_id])));

      first = last;
    }

    std::for_each(collect_futures.begin(), collect_futures.end(),
                  std::mem_fn(&std::future<void>::wait));
  }

  return dst;
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

[[nodiscard]] static auto FindWindows(
    State& state, std::vector<ReadOverlapsPair> const& reads_overlaps,
    std::unique_ptr<biosoup::NucleicAcid> const& query_read,
    double const match_ratio) -> std::vector<PloidyInterval> {
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

[[nodiscard]] static auto CorrectReads(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>::iterator seq_first,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>::iterator seq_last,
    std::vector<PloidyInterval>::iterator interval_firsts,
    std::vector<PloidyInterval>::iterator interval_last)
    -> std::vector<AnnotatedRead> {
  auto dst = std::vector<AnnotatedRead>();

  return dst;
}

}  // namespace detail

auto SnpErrorCorrect(
    State& state, std::vector<std::unique_ptr<biosoup::NucleicAcid>> src_reads)
    -> std::vector<AnnotatedRead> {
  auto dst = std::vector<AnnotatedRead>();

  auto timer = biosoup::Timer();

  timer.Start();
  auto reads_overlaps =
      detail::FindOverlapsAndFilterReads(state, std::move(src_reads));
  fmt::print(stderr,
             "[camel::SnpErrorCorrect]({:12.3f}) tied reads with overlaps\n",
             timer.Stop());

  auto window_futures = std::deque<
      std::future<std::vector<std::vector<detail::PloidyInterval>>>>();

  auto const can_pop_future =
      [](std::future<std::vector<std::vector<detail::PloidyInterval>>>& f)
      -> bool {
    return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
  };

  {
    auto cnt = 0U;
    timer.Start();
    for (auto first = reads_overlaps.cbegin();
         first != reads_overlaps.cend();) {
      auto const last =
          std::min(std::next(first, detail::kPackSize), reads_overlaps.cend());

      window_futures.emplace_back(state.thread_pool->Submit(
          [&state, &reads_overlaps](
              std::vector<ReadOverlapsPair>::const_iterator first,
              std::vector<ReadOverlapsPair>::const_iterator last) {
            auto dst = std::vector<std::vector<detail::PloidyInterval>>();
            dst.reserve(std::distance(first, last));

            std::transform(first, last, std::back_inserter(dst),
                           [&state, &reads_overlaps](ReadOverlapsPair const& ro)
                               -> std::vector<detail::PloidyInterval> {
                             return detail::FindWindows(state, reads_overlaps,
                                                        ro.read, 0.9);
                           });

            return dst;
          },
          first, last));

      first = last;
    }

    while (!window_futures.empty()) {
      auto const windows = window_futures.front().get();
      window_futures.pop_front();

      cnt += windows.size();
      if (cnt % 1000 == 0U || cnt == reads_overlaps.size()) {
        fmt::print(stderr, "\r[camel::SnpErrorCorrect]({:12.3f}) {} / {}",
                   timer.Lap(), cnt, reads_overlaps.size());
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
