#include "camel/mapping.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <deque>
#include <iterator>
#include <numeric>

#include "biosoup/timer.hpp"
#include "detail/overlap.h"
#include "fmt/core.h"
#include "ram/minimizer_engine.hpp"
#include "tsl/robin_set.h"

namespace camel {

namespace detail {

static constexpr auto kMinimizeBatchCap = 1UL << 34UL;
static constexpr auto kMapBatchCap = 1UL << 26UL;

static constexpr auto kDefaultExpectedCoverage = 32U;
static constexpr auto kDefaultWindowLen = 500U;

static auto FindBatchLast(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator first,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator last,
    std::uint64_t batch_cap)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator {
  for (auto batch_sz = 0UL; batch_sz < batch_cap && first != last; ++first) {
    batch_sz += (*first)->inflated_len;
  }

  return first;
}

static auto SweepOverlaps(std::vector<biosoup::Overlap>& overlaps,
                          std::uint32_t const kReadLen,
                          std::uint32_t const kWindowLen,
                          std::uint32_t const kExpectedCoverage)
    -> std::uint64_t {
  auto windows =
      std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>>();
  windows.resize(std::round(1.0 * kReadLen / kWindowLen));

  for (auto& window : windows) {
    windows.reserve(kExpectedCoverage + 1U);
  }

  for (auto ovlp_id = 0U; ovlp_id < overlaps.size(); ++ovlp_id) {
    auto first_idx = overlaps[ovlp_id].lhs_begin / kWindowLen;
    auto last_idx =
        std::min(static_cast<std::uint64_t>(
                     std::ceil(1.0 * overlaps[ovlp_id].lhs_end / kWindowLen)) +
                     1U,
                 windows.size());

    for (; first_idx != last_idx; ++first_idx) {
      windows[first_idx].emplace_back(overlaps[ovlp_id].score, ovlp_id);
      for (auto i = windows[first_idx].size() - 1U; i > 0U; --i) {
        if (windows[first_idx][i - 1].first < windows[first_idx][i].first) {
          std::swap(windows[first_idx][i - 1], windows[first_idx][i]);
        }
      }

      while (windows[first_idx].size() > kExpectedCoverage) {
        windows[first_idx].resize(kExpectedCoverage);
      }
    }
  }

  auto valid_ids = tsl::robin_set<std::uint32_t>();
  valid_ids.reserve(windows.size() * kExpectedCoverage);

  for (auto const& window : windows) {
    for (auto const [score, idx] : window) {
      valid_ids.insert(idx);
    }
  }

  auto updated_ovlps = std::vector<biosoup::Overlap>();
  updated_ovlps.reserve(valid_ids.size());

  for (auto const idx : valid_ids) {
    updated_ovlps.push_back(overlaps[idx]);
  }

  auto dst = overlaps.size() - updated_ovlps.size();
  std::swap(overlaps, updated_ovlps);

  return dst;
}

auto FormatAdjacencyList(
    State& state,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads,
    std::vector<std::vector<biosoup::Overlap>> overlaps)
    -> std::vector<std::vector<biosoup::Overlap>> {
  auto ovlp_futures = std::vector<std::future<void>>();
  ovlp_futures.reserve(reads.size());

  auto const transform_target =
      [&overlaps](std::uint32_t const target_id) -> void {
    auto const target_cnt =
        std::count_if(overlaps[target_id].cbegin(), overlaps[target_id].cend(),
                      [target_id](biosoup::Overlap const& ovlp) -> bool {
                        return ovlp.rhs_id < target_id;
                      });

    auto dst = std::vector<biosoup::Overlap>();
    dst.reserve(target_cnt);

    for (auto const& ovlp : overlaps[target_id]) {
      if (ovlp.rhs_id < target_id) {
        dst.push_back(camel::detail::ReverseOverlap(ovlp));
      }
    }

    auto const by_query_id_cmp = [](biosoup::Overlap const& a,
                                    biosoup::Overlap const& b) -> bool {
      return a.lhs_id < b.lhs_id;
    };

    std::sort(dst.begin(), dst.end(), by_query_id_cmp);
    std::swap(overlaps[target_id], dst);
  };

  auto const transform_range = [&transform_target](
                                   std::uint32_t const first_id,
                                   std::uint32_t const last_id) -> void {
    for (auto target_id = first_id; target_id != last_id; ++target_id) {
      transform_target(target_id);
    }
  };

  for (auto target_id = 0U; target_id < reads.size(); target_id += 1000U) {
    ovlp_futures.emplace_back(state.thread_pool->Submit(
        transform_range, target_id,
        std::min(target_id + 1000U, static_cast<std::uint32_t>(reads.size()))));
  }

  std::for_each(ovlp_futures.begin(), ovlp_futures.end(),
                std::mem_fn(&std::future<void>::get));

  return overlaps;
}

}  // namespace detail

MapCfg::MapCfg(std::uint8_t kmer_len, std::uint8_t win_len, double filter_p)
    : kmer_len(kmer_len), win_len(win_len), filter_p(filter_p) {}

auto FindOverlaps(
    State& state, MapCfg const map_cfg,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads)
    -> std::vector<std::vector<biosoup::Overlap>> {
  for (auto idx = 1U; idx < reads.size(); ++idx) {
    if (reads[idx - 1U]->id + 1U != reads[idx]->id) {
      throw std::runtime_error(
          "[camel::FindOverlaps] read ids must form continuous ascending "
          "sequence");
    }
  }

  auto read_ovlps = std::vector<std::vector<biosoup::Overlap>>(reads.size());

  auto minimizer_engine = ram::MinimizerEngine(
      state.thread_pool, map_cfg.kmer_len, map_cfg.win_len);

  auto map_futures = std::vector<std::future<std::vector<biosoup::Overlap>>>();

  auto const store_ovlp = [&read_ovlps](biosoup::Overlap const& ovlp) -> void {
    auto& target_vec = read_ovlps[ovlp.rhs_id];
    if (target_vec.size() == target_vec.capacity()) {
      target_vec.reserve(target_vec.size() * 1.2);
    }

    target_vec.push_back(ovlp);
  };

  auto const find_batch_batckwards =
      [](std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator
             first,
         std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator
             last,
         std::size_t const batch_cap)
      -> std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator {
    if (last != first) {
      auto curr_read = std::prev(last);
      for (auto batch_sz = 0UL; batch_sz < batch_cap;
           std::advance(curr_read, -1)) {
        batch_sz += (*curr_read)->inflated_len;
        if (curr_read == first) {
          break;
        }
      }

      return curr_read;
    } else {
      return first;
    }
  };

  auto const find_batch_forward =
      [](std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator
             first,
         std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator
             last,
         std::size_t const batch_cap)
      -> std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator {
    for (auto batch_sz = 0UL; batch_sz < batch_cap && first != last; ++first) {
      batch_sz += (*first)->inflated_len;
    }

    return first;
  };

  auto const map_sequence =
      [&minimizer_engine](std::unique_ptr<biosoup::NucleicAcid> const& read)
      -> std::vector<biosoup::Overlap> {
    return minimizer_engine.Map(read, true, true, true);
  };

  auto const map_sequence_async =
      [&thread_pool = state.thread_pool, &map_sequence](
          std::reference_wrapper<std::unique_ptr<biosoup::NucleicAcid> const>
              read) -> std::future<std::vector<biosoup::Overlap>> {
    return thread_pool->Submit(map_sequence, read);
  };

  auto defrag_futures = std::deque<std::future<void>>();

  auto const is_defrag_ready = [](std::future<void> const& df) -> bool {
    return df.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
  };

  auto timer = biosoup::Timer();
  auto minimize_last = reads.cend();
  while (true) {
    timer.Start();

    while (!defrag_futures.empty() && is_defrag_ready(defrag_futures.front())) {
      defrag_futures.front().wait();
      defrag_futures.pop_front();
    }

    auto const minimize_first = find_batch_batckwards(
        reads.cbegin(), minimize_last, detail::kMinimizeBatchCap);

    minimizer_engine.Minimize(minimize_first, minimize_last, true);
    minimizer_engine.Filter(map_cfg.filter_p);

    fmt::print(
        stderr, "[camel::FindOverlaps]({:12.3f}) minimized {} / {} reads\n",
        timer.Stop(), std::distance(minimize_first, reads.end()), reads.size());

    timer.Start();
    for (auto map_first = reads.cbegin(); map_first < minimize_last;) {
      auto const map_last =
          find_batch_forward(map_first, minimize_last, detail::kMapBatchCap);

      map_futures.reserve(std::distance(map_first, map_last));
      std::transform(map_first, map_last, std::back_inserter(map_futures),
                     map_sequence_async);

      for (auto& it : map_futures) {
        for (auto&& ovlp : it.get()) {
          store_ovlp(ovlp);
        }
      }

      fmt::print(stderr,
                 "\r[camel::FindOverlaps]({:12.3f}) mapped {} / {} reads",
                 timer.Lap(), std::distance(reads.cbegin(), map_last),
                 std::distance(reads.begin(), minimize_last));

      map_first = map_last;
      map_futures.clear();
    }

    timer.Stop();
    fmt::print(stderr, "\n");

    defrag_futures.emplace_back(state.thread_pool->Submit(
        [&read_ovlps](
            std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator
                first,
            std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator
                last) -> void {
          for (auto it = first; it != last; ++it) {
            auto& ovlp_vec = read_ovlps[(*it)->id];
            if (ovlp_vec.size() < ovlp_vec.capacity()) {
              std::remove_reference_t<decltype(ovlp_vec)>(ovlp_vec.cbegin(),
                                                          ovlp_vec.cend())
                  .swap(ovlp_vec);
            }
          }
        },
        minimize_first, minimize_last));

    if (minimize_first != reads.cbegin()) {
      minimize_last = std::prev(minimize_first);
    } else {
      break;
    }
  }

  timer.Start();
  while (!defrag_futures.empty()) {
    defrag_futures.front().wait();
    defrag_futures.pop_front();
  }
  timer.Stop();

  return read_ovlps;
}

auto FindConfidentOverlaps(
    State& state, MapCfg const map_cfg,
    std::uint32_t const kWindowLen,
    std::uint32_t const kMaxCoverage,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& src_reads)
    -> std::vector<std::vector<biosoup::Overlap>> {
  auto overlaps = std::vector<std::vector<biosoup::Overlap>>(src_reads.size());
  auto sweep_futures = std::vector<std::future<std::uint64_t>>();

  auto minimizer_engine = ram::MinimizerEngine(
      state.thread_pool, map_cfg.kmer_len, map_cfg.win_len);
  auto map_futures = std::vector<std::future<std::vector<biosoup::Overlap>>>();

  auto const store_ovlp = [&src_reads,
                           &overlaps](biosoup::Overlap const& ovlp) -> void {
    if (detail::DetermineOverlapType(ovlp, src_reads[ovlp.lhs_id]->inflated_len,
                                     src_reads[ovlp.rhs_id]->inflated_len) >
        detail::OverlapType::kUnclassified) {
      overlaps[ovlp.lhs_id].push_back(ovlp);
      overlaps[ovlp.rhs_id].push_back(detail::ReverseOverlap(ovlp));
    }
  };

  auto const map_sequence =
      [&minimizer_engine](std::unique_ptr<biosoup::NucleicAcid> const& read)
      -> std::vector<biosoup::Overlap> {
    return minimizer_engine.Map(read, true, true, true);
  };

  auto const map_sequence_async =
      [&thread_pool = state.thread_pool, &map_sequence](
          std::reference_wrapper<std::unique_ptr<biosoup::NucleicAcid> const>
              read) -> std::future<std::vector<biosoup::Overlap>> {
    return thread_pool->Submit(map_sequence, read);
  };

  auto timer = biosoup::Timer();
  for (auto minimize_first = src_reads.begin();
       minimize_first != src_reads.end();) {
    timer.Start();
    auto const minimize_last = detail::FindBatchLast(
        minimize_first, src_reads.end(), detail::kMinimizeBatchCap);

    minimizer_engine.Minimize(minimize_first, minimize_last, true);
    minimizer_engine.Filter(map_cfg.filter_p);

    fmt::print(
        stderr,
        "[camel::FindConfidentOverlaps]({:12.3f}) minimized {} / {} reads\n",
        timer.Stop(), std::distance(src_reads.begin(), minimize_last),
        src_reads.size());

    timer.Start();
    for (auto map_first = src_reads.begin(); map_first != minimize_last;) {
      auto const map_last =
          detail::FindBatchLast(map_first, minimize_last, detail::kMapBatchCap);

      map_futures.reserve(std::distance(map_first, map_last));
      std::transform(map_first, map_last, std::back_inserter(map_futures),
                     map_sequence_async);

      for (auto& it : map_futures) {
        for (auto const& ovlp : it.get()) {
          store_ovlp(ovlp);
        }
      }

      map_futures.clear();

      {
        auto first = (*minimize_first)->id < (*map_first)->id ? minimize_first
                                                              : map_first;
        auto last = minimize_last;

        for (; first != last; ++first) {
          if (overlaps[(*first)->id].size() > 420 * 4) { // NOTE: this might need re adjusting
            sweep_futures.emplace_back(state.thread_pool->Submit(
                [&src_reads,
                 &overlaps](std::uint32_t const id) -> std::uint64_t {
                  return detail::SweepOverlaps(
                      overlaps[id], src_reads[id]->inflated_len, 320, 40);
                },
                (*first)->id));
          }
        }

        auto const n_discards = std::transform_reduce(
            sweep_futures.begin(), sweep_futures.end(), 0U,
            std::plus<std::uint64_t>(),
            std::mem_fn(&std::future<std::uint64_t>::get));

        sweep_futures.clear();

        // fmt::print(
        //     stderr,
        //     "\r[camel::FindConfidentOverlaps]({:12.3f}) swept {} low quality
        //     " "overlaps", timer.Stop(), n_discards);
      }

      fmt::print(
          stderr,
          "\r[camel::FindConfidentOverlaps]({:12.3f}) mapped {} / {} reads",
          timer.Lap(), std::distance(src_reads.begin(), map_last),
          std::distance(src_reads.begin(), minimize_last));

      map_first = map_last;
    }

    timer.Stop();
    fmt::print(stderr, "\n");
    minimize_first = minimize_last;
  }

  timer.Start();
  overlaps = detail::FormatAdjacencyList(state, src_reads, std::move(overlaps));
  fmt::print(
      stderr,
      "[camel::FindConfidentOverlaps]({:12.3f}) formated adjacency list\n",
      timer.Stop());

  return overlaps;
}

}  // namespace camel
