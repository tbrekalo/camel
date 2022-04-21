#include "camel/mapping.h"

#include <algorithm>
#include <chrono>
#include <deque>
#include <iterator>
#include <numeric>

#include "biosoup/progress_bar.hpp"
#include "biosoup/timer.hpp"
#include "detail/overlap.h"
#include "fmt/core.h"
#include "fmt/ostream.h"
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
    -> std::vector<ReadOverlapIndex> {
  auto read_ovlps = std::vector<ReadOverlapIndex>(reads.size());

  auto minimizer_engine =
      ram::MinimizerEngine(thread_pool, map_cfg.kmer_len, map_cfg.win_len);

  auto map_futures = std::vector<std::future<std::vector<biosoup::Overlap>>>();

  auto const store_ovlp = [&read_ovlps](biosoup::Overlap const& ovlp) -> void {
    auto& target_vec = read_ovlps[ovlp.rhs_id].target_overlaps;
    if (target_vec.size() == target_vec.capacity()) {
      target_vec.reserve(target_vec.size() * 1.5);
    }

    target_vec.push_back(ovlp);
  };

  auto const find_batch_last =
      [](std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator
             first,
         std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator
             last,
         std::size_t const batch_cap) {
        for (auto batch_sz = 0UL; batch_sz < batch_cap && first != last;
             ++first) {
          batch_sz += first->get()->inflated_len;
        }

        return first;
      };

  auto const map_sequence =
      [&minimizer_engine](std::unique_ptr<biosoup::NucleicAcid> const& read)
      -> std::vector<biosoup::Overlap> {
    auto dst = minimizer_engine.Map(read, true, true, true);
    std::sort(  // TODO: does ram provide sort order?
        dst.begin(), dst.end(),
        [](biosoup::Overlap const& a, biosoup::Overlap const& b) -> bool {
          return a.rhs_id < b.rhs_id;
        });

    return dst;
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
  for (auto minimize_first = reads.cbegin(); minimize_first != reads.cend();) {
    timer.Start();

    while (!defrag_futures.empty() && is_defrag_ready(defrag_futures.front())) {
      defrag_futures.front().wait();
      defrag_futures.pop_front();
    }

    auto const minimize_last = find_batch_last(minimize_first, reads.cend(),
                                               detail::kMinimizeBatchCap);

    minimizer_engine.Minimize(minimize_first, minimize_last, true);
    minimizer_engine.Filter(map_cfg.filter_p);

    fmt::print(stderr,
               "[camel::FindOverlaps]({:12.3f}) minimized {} / {} reads\n",
               timer.Stop(), std::distance(reads.cbegin(), minimize_last),
               reads.size());

    for (auto map_first = minimize_first; map_first < reads.cend();) {
      timer.Start();
      auto const map_last =
          find_batch_last(map_first, reads.cend(), detail::kMapBatchCap);

      map_futures.reserve(std::distance(map_first, map_last));
      std::transform(map_first, map_last, std::back_inserter(map_futures),
                     map_sequence_async);

      for (auto& it : map_futures) {
        auto ovlp_vec = it.get();
        if (!ovlp_vec.empty()) {
          auto const query_id = ovlp_vec.front().lhs_id;
          auto& query_remotes = read_ovlps[query_id].remote_intervals;

          for (auto target_first = ovlp_vec.cbegin();
               target_first != ovlp_vec.cend();) {
            auto const target_id = target_first->rhs_id;
            auto const target_last = std::find_if_not(
                target_first, ovlp_vec.cend(),
                [target_id](biosoup::Overlap const& ovlp) -> bool {
                  return ovlp.rhs_id == target_id;
                });

            auto const& target_vec = read_ovlps[target_id].target_overlaps;

            query_remotes.push_back(OverlapInterval{
                .target_id = target_id, .first_index = target_vec.size()});
            std::for_each(target_first, target_last, store_ovlp);
            query_remotes.back().last_index = target_vec.size();

            target_first = target_last;
          }
        }
      }

      fmt::print(stderr,
                 "[camel::FindOverlaps]({:12.3f}) mapped {} / {} reads\n",
                 timer.Stop(), std::distance(reads.cbegin(), map_last),
                 std::distance(reads.cbegin(), minimize_last));

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
            read_ovlps[(*it)->id].target_overlaps.shrink_to_fit();
            read_ovlps[(*it)->id].remote_intervals.shrink_to_fit();
          }
        },
        minimize_first, minimize_last));

    minimize_first = minimize_last;
  }

  timer.Start();
  while (!defrag_futures.empty()) {
    defrag_futures.front().wait();
    defrag_futures.pop_front();
  }
  timer.Stop();

  fmt::print(stderr, "[camel::FindOverlaps]({:12.3f})\n", timer.elapsed_time());
  return read_ovlps;
}

}  // namespace camel
