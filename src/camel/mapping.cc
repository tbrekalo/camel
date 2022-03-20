#include "camel/mapping.h"

#include <algorithm>
#include <iterator>

#include "biosoup/timer.hpp"
#include "detail/overlap.h"
#include "fmt/core.h"
#include "ram/minimizer_engine.hpp"

namespace camel {

static constexpr std::size_t kMinimizeBatchCap = 1UL << 32UL;
static constexpr std::size_t kMapBatchCap = 1UL << 30UL;

auto FindOverlaps(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& seqs,
    MapCfg map_cfg) -> std::vector<std::vector<biosoup::Overlap>> {
  auto dst = std::vector<std::vector<biosoup::Overlap>>(
      seqs.size(), std::vector<biosoup::Overlap>());

  auto store_overlap = [&dst](biosoup::Overlap const& ovlp) -> void {
    dst[ovlp.lhs_id].push_back(ovlp);
    dst[ovlp.rhs_id].push_back(detail::ReverseOverlap(ovlp));
  };

  auto minimizer_engine =
      ram::MinimizerEngine(thread_pool, map_cfg.kmer_len, map_cfg.win_len);

  auto map_futures = std::vector<std::future<std::vector<biosoup::Overlap>>>();

  auto const map_sequence =
      [&minimizer_engine](std::unique_ptr<biosoup::NucleicAcid> const& seq)
      -> std::vector<biosoup::Overlap> {
    return minimizer_engine.Map(seq, true, true);
  };

  auto const async_map_sequence =
      [&thread_pool, &map_sequence](
          std::reference_wrapper<std::unique_ptr<biosoup::NucleicAcid> const>
              seq) { return thread_pool->Submit(map_sequence, seq); };

  auto const find_batch_end =
      [&seqs](std::size_t const begin_idx, std::size_t const end_idx,
              std::size_t const batch_cap) -> std::size_t {
    auto curr_idx = begin_idx;
    for (auto batch_sz = 0UL; curr_idx < end_idx && batch_sz < batch_cap;
         ++curr_idx) {
      batch_sz += seqs[curr_idx]->inflated_len;
    }

    return curr_idx;
  };

  auto timer = biosoup::Timer();
  timer.Start();

  for (auto minimize_batch_begin = 0UL; minimize_batch_begin < seqs.size();) {
    auto const minimize_batch_end =
        find_batch_end(minimize_batch_begin, seqs.size(), kMinimizeBatchCap);

    minimizer_engine.Minimize(std::next(seqs.begin(), minimize_batch_begin),
                              std::next(seqs.begin(), minimize_batch_end));

    minimizer_engine.Filter(map_cfg.filter_p);

    fmt::print(stderr,
               "[camel::FindOverlaps]({:12.3f}) minimized {} sequences\n",
               timer.Stop(), minimize_batch_end - minimize_batch_begin);
    timer.Start();

    for (auto map_batch_begin = 0U; map_batch_begin < minimize_batch_end;) {
      auto const map_batch_end =
          find_batch_end(map_batch_begin, minimize_batch_end, kMapBatchCap);

      map_futures.reserve(map_batch_end - map_batch_begin);

      std::transform(std::next(seqs.cbegin(), map_batch_begin),
                     std::next(seqs.cbegin(), map_batch_end),
                     std::back_inserter(map_futures), async_map_sequence);

      for (auto& it : map_futures) {
        for (auto&& ovlp : it.get()) {
          store_overlap(ovlp);
        }
      }

      fmt::print(stderr,
                 "[camel::FindOverlaps]({:12.3f}) mapped {} sequences\n",
                 timer.Stop(), map_batch_end - map_batch_begin);
      timer.Start();

      map_batch_begin = map_batch_end;
      map_futures.clear();
    }

    minimize_batch_begin = minimize_batch_end;
  }
  
  auto const n_ovlps = std::accumulate(
      dst.cbegin(), dst.cend(), 0UL,
      [](std::size_t const init, std::vector<biosoup::Overlap> const& vec)
          -> std::size_t { 
            return init + vec.size(); 
          });

  fmt::print(stderr, "[camel::FindOverlaps] found {} overlaps\n", n_ovlps);
  return dst;
}

}  // namespace camel
