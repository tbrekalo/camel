#include "camel/mapping.h"

#include <algorithm>
#include <chrono>
#include <deque>
#include <iterator>
#include <numeric>

#include "biosoup/timer.hpp"
#include "detail/overlap.h"
#include "fmt/core.h"
#include "ram/minimizer_engine.hpp"

namespace camel {

namespace detail {

static constexpr auto kMinimizeBatchCap = 1UL << 32UL;
static constexpr auto kMapBatchCap = 1UL << 30UL;

static constexpr auto kCoverageFactor = 12UL;

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

struct SweepEvent {
  std::uint32_t pos;
  std::uint32_t ovlp_id;
};

static auto SweepOverlaps(std::vector<biosoup::Overlap>& overlaps,
                          std::uint64_t estimated_covg)
    -> std::pair<std::uint64_t, std::uint64_t> {
  auto const cmp_lambda =
      [](std::pair<std::uint32_t, std::uint32_t> const& lhs_pis,
         std::pair<std::uint32_t, std::uint32_t> const& rhs_pis) -> bool {
    return lhs_pis.second < rhs_pis.second;
  };

  auto active_sz = 0U;
  auto active_buff = std::array<std::pair<std::uint32_t, std::uint32_t>,
                                detail::kCoverageFactor>();

  auto events = std::vector<SweepEvent>();

  events.reserve(overlaps.size() * 2U);
  for (auto i = 0U; i < overlaps.size(); ++i) {
    events.emplace_back(SweepEvent{.pos = overlaps[i].lhs_begin, .ovlp_id = i});
    events.emplace_back(SweepEvent{.pos = overlaps[i].lhs_end, .ovlp_id = i});
  }

  std::sort(events.begin(), events.end(),
            [](SweepEvent const& lhs, SweepEvent const& rhs) -> bool {
              return lhs.pos < rhs.pos;
            });

  for (auto const& event : events) {
    if (active_sz > 0U) {
      auto const active_end = std::next(active_buff.begin(), active_sz);
      auto opt_iter = std::find_if(
          active_buff.begin(), active_end,
          [qo_id = event.ovlp_id](std::pair<std::uint32_t, std::uint32_t> pis)
              -> bool { return pis.first == qo_id; });

      if (opt_iter != active_buff.end()) {
        auto last_iter = std::prev(active_end);
        if (opt_iter != last_iter) {
          std::swap(*opt_iter, *last_iter);
        }

        --active_sz;
      }
    }

    if (active_sz < active_buff.size()) {
      active_buff[active_sz++] =
          std::pair(event.ovlp_id, overlaps[event.ovlp_id].score);
    } else {
      auto min_iter = std::min_element(
          active_buff.begin(), std::next(active_buff.begin(), active_sz),
          cmp_lambda);

      if (overlaps[min_iter->first].score < overlaps[event.ovlp_id].score) {
        *min_iter = std::pair(event.ovlp_id, overlaps[event.ovlp_id].score);
      } else {
        std::swap(overlaps[event.ovlp_id].lhs_begin,
                  overlaps[event.ovlp_id].lhs_end);
      }
    }
  }

  auto remove_iter = std::remove_if(overlaps.begin(), overlaps.end(),
                                    [](biosoup::Overlap const& ovlp) -> bool {
                                      return ovlp.lhs_begin > ovlp.rhs_end;
                                    });

  auto dst_cnt = std::distance(remove_iter, overlaps.end());
  for (auto curr = remove_iter; curr != overlaps.end(); ++curr) {
    estimated_covg -= curr->lhs_begin - curr->lhs_end;
  }

  overlaps.erase(remove_iter, overlaps.cend());
  return {estimated_covg, dst_cnt};
}

}  // namespace detail

MapCfg::MapCfg(std::uint8_t kmer_len, std::uint8_t win_len, double filter_p)
    : kmer_len(kmer_len), win_len(win_len), filter_p(filter_p) {}

auto FindOverlaps(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool, MapCfg const map_cfg,
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

  auto minimizer_engine =
      ram::MinimizerEngine(thread_pool, map_cfg.kmer_len, map_cfg.win_len);

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
      [&thread_pool, &map_sequence](
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

    for (auto map_first = reads.cbegin(); map_first < minimize_last;) {
      timer.Start();
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
                 "[camel::FindOverlaps]({:12.3f}) mapped {} / {} reads\n",
                 timer.Stop(), std::distance(reads.cbegin(), map_last),
                 std::distance(reads.begin(), minimize_last));

      map_first = map_last;
      map_futures.clear();
    }

    defrag_futures.emplace_back(thread_pool->Submit(
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
    std::shared_ptr<thread_pool::ThreadPool> thread_pool, MapCfg const map_cfg,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& src_reads)
    -> std::vector<std::vector<biosoup::Overlap>> {
  auto overlaps = std::vector<std::vector<biosoup::Overlap>>(src_reads.size());
  auto estimated_covgs = std::vector<std::uint64_t>(src_reads.size());

  auto sweep_futures = std::vector<std::future<std::uint64_t>>();

  auto minimizer_engine =
      ram::MinimizerEngine(thread_pool, map_cfg.kmer_len, map_cfg.win_len);
  auto map_futures = std::vector<std::future<std::vector<biosoup::Overlap>>>();

  auto const store_ovlp = [src_reads, &overlaps, &estimated_covgs](
                              biosoup::Overlap const& ovlp) -> void {
    if (detail::DetermineOverlapType(ovlp, src_reads[ovlp.lhs_id]->inflated_len,
                                     src_reads[ovlp.rhs_id]->inflated_len) >
        detail::OverlapType::kUnclassified) {
      estimated_covgs[ovlp.lhs_id] += ovlp.lhs_end - ovlp.lhs_end;
      estimated_covgs[ovlp.rhs_id] += ovlp.rhs_end - ovlp.rhs_begin;

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
      [&thread_pool, &map_sequence](
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
    minimizer_engine.Filter(0.001);

    fmt::print(stderr,
               "[camel::SnpErrorCorrect]({:12.3f}) minimized {} / {} reads\n",
               timer.Stop(), std::distance(src_reads.begin(), minimize_last),
               src_reads.size());

    for (auto map_first = src_reads.begin(); map_first != minimize_last;) {
      timer.Start();
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

      fmt::print(stderr,
                 "[camel::SnpErrorCorrect]({:12.3f}) mapped {} / {} reads\n",
                 timer.Stop(), std::distance(src_reads.begin(), map_last),
                 std::distance(src_reads.begin(), minimize_last));

      {
        timer.Start();
        auto first = (*minimize_first)->id < (*map_first)->id ? minimize_first
                                                              : map_first;
        auto last = minimize_last;

        for (; first != last; ++first) {
          if (estimated_covgs[(*first)->id] >=
              detail::kCoverageFactor * (*first)->inflated_len) {
            sweep_futures.emplace_back(thread_pool->Submit(
                [&overlaps,
                 &estimated_covgs](std::uint32_t const id) -> std::uint64_t {
                  auto const [estimated_covg, discard_cnt] =
                      detail::SweepOverlaps(overlaps[id], estimated_covgs[id]);

                  estimated_covgs[id] = estimated_covg;
                  return discard_cnt;
                },
                (*first)->id));
          }
        }

        auto const n_discards = std::transform_reduce(
            sweep_futures.begin(), sweep_futures.end(), 0U,
            std::plus<std::uint64_t>(),
            std::mem_fn(&std::future<std::uint64_t>::get));

        sweep_futures.clear();

        fmt::print(stderr,
                   "[camel::SnpErrorCorrect]({:12.3f}) swept {} low quality "
                   "overlaps\n",
                   timer.Stop(), n_discards);
      }

      map_first = map_last;
    }

    minimize_first = minimize_last;
  }

  return overlaps;
}

}  // namespace camel
