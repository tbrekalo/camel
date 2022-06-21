#include "camel/mapping.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <deque>
#include <iterator>
#include <numeric>
#include <optional>

#include "biosoup/timer.hpp"
#include "detail/overlap.h"
#include "fmt/core.h"
#include "ram/minimizer_engine.hpp"
#include "tsl/robin_set.h"

namespace camel {

namespace detail {

static constexpr auto kMinimizeBatchCap = 1UL << 34UL;
static constexpr auto kMapBatchCap = 1UL << 26UL;

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

  auto const find_batch_forwards =
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

    auto const minimize_first = find_batch_forwards(
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
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads)
    -> std::vector<std::vector<biosoup::Overlap>> {
  auto dst = std::vector<std::vector<biosoup::Overlap>>(reads.size());
  auto primaries = std::vector<biosoup::Overlap>(reads.size());

  auto minimizer_engine = ram::MinimizerEngine(
      state.thread_pool, map_cfg.kmer_len, map_cfg.win_len);

  auto const ovlp_len = [](biosoup::Overlap const& ovlp) -> std::uint32_t {
    return std::max(ovlp.lhs_end - ovlp.lhs_begin,
                    ovlp.rhs_end - ovlp.lhs_begin);
  };

  auto const ovlp_err = [&ovlp_len](biosoup::Overlap const& ovlp) -> double {
    return 1.0 - static_cast<double>(std::min(ovlp.lhs_end - ovlp.lhs_begin,
                                              ovlp.rhs_end - ovlp.rhs_begin)) /
                     static_cast<double>(ovlp_len(ovlp));
  };

  auto const find_batch_last =
      [](std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator
             first,
         std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator
             last,
         std::uint64_t const kBatchCap)
      -> std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator {
    for (auto batch_sz = 0UL; batch_sz < kBatchCap && first != last; ++first) {
      batch_sz += (*first)->inflated_len;
    }

    return first;
  };

  auto timer = biosoup::Timer();

  {
    auto map_futures =
        std::vector<std::future<std::optional<biosoup::Overlap>>>();
    for (auto minimize_first = reads.cbegin();
         minimize_first != reads.cend();) {
      timer.Start();
      auto const minimize_last = find_batch_last(minimize_first, reads.cend(),
                                                 detail::kMinimizeBatchCap);

      minimizer_engine.Minimize(minimize_first, minimize_last, true);
      fmt::print(stderr,
                 "[camel::FindConfidentOverlaps]({:12.3f}) minimized "
                 "{} / {} reads\n",
                 timer.Stop(), std::distance(reads.cbegin(), minimize_last),
                 reads.size());

      timer.Start();
      for (auto map_first = minimize_first; map_first != minimize_last;) {
        auto const map_last =
            find_batch_last(map_first, minimize_last, detail::kMapBatchCap);

        map_futures.reserve(std::distance(map_first, map_last));
        for (auto iter = map_first; iter != map_last; ++iter) {
          map_futures.emplace_back(state.thread_pool->Submit(
              [&minimizer_engine, &ovlp_len,
               &ovlp_err](std::vector<
                          std::unique_ptr<biosoup::NucleicAcid>>::const_iterator
                              read_iter) -> std::optional<biosoup::Overlap> {
                auto ovlps = minimizer_engine.Map(*read_iter, true, true);
                ovlps.erase(
                    std::remove_if(
                        ovlps.begin(), ovlps.end(),
                        [&ovlp_err](biosoup::Overlap const& ovlp) -> bool {
                          return ovlp_err(ovlp) > 0.2;
                        }),

                    ovlps.end());

                if (!ovlps.empty()) {
                  return *std::max_element(
                      ovlps.cbegin(), ovlps.cend(),
                      [&ovlp_len](biosoup::Overlap const& a,
                                  biosoup::Overlap const& b) -> bool {
                        return ovlp_len(a) < ovlp_len(b);
                      });
                } else {
                  return std::nullopt;
                }
              },
              iter));
        }

        for (auto& future : map_futures) {
          if (auto const opt_ovlp = future.get(); opt_ovlp) {
            if (ovlp_len(primaries[opt_ovlp->lhs_id]) < ovlp_len(*opt_ovlp)) {
              primaries[opt_ovlp->lhs_id] = *opt_ovlp;
            }
          }
        }

        fmt::print(
            stderr,
            "\r[camel::FindConfidentOverlaps]({:12.3f}) mapped {} / {} reads",
            timer.Lap(), std::distance(minimize_first, map_last),
            std::distance(minimize_first, minimize_last));

        map_futures.clear();
        map_first = map_last;
      }

      timer.Stop();
      fmt::print(stderr, "\n");

      minimize_first = minimize_last;
    }
  }

  {
    timer.Start();

    for (auto const& ovlp : primaries) {
      dst[ovlp.rhs_id].push_back(detail::ReverseOverlap(ovlp));
    }

    auto target_cnt = 0U;
    auto query_cnt = 0U;

    for (auto& ovlps : dst) {
      if (!ovlps.empty()) {
        std::vector<biosoup::Overlap>(ovlps).swap(ovlps);
        ++target_cnt;
      } else {
        ++query_cnt;
      }
    }

    fmt::print(stderr,
               "[camel:FindConfidentOverlaps]({:12.3f}) organized overlaps in "
               "{} targets and {} queries\n",
               timer.Stop(), target_cnt, query_cnt);
  }

  return dst;
}

}  // namespace camel
