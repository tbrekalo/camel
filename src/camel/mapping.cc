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

}  // namespace camel
