#include "camel/mapping.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <numeric>

#include "biosoup/timer.hpp"
#include "detail/overlap.h"
#include "fmt/core.h"
#include "ram/minimizer_engine.hpp"

namespace camel {

namespace detail {

static constexpr std::size_t kMapBatchCap = 1UL << 30UL;

}

MapCfg::MapCfg(std::uint8_t kmer_len, std::uint8_t win_len, double filter_p)
    : kmer_len(kmer_len), win_len(win_len), filter_p(filter_p) {}

auto FindOverlaps(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool, MapCfg const map_cfg,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& seqs,
    std::size_t const target_group_sz) -> std::vector<OverlapGroup> {
  auto dst = std::vector<OverlapGroup>();
  auto seq_ovlp_group_ids = std::vector<std::uint32_t>(seqs.size(), 0U);

  auto const find_batch_last =
      [&seqs](std::uint32_t const first_idx, std::uint32_t const last_idx,
              std::size_t const batch_cap) -> std::uint32_t {
    auto curr_idx = first_idx;
    for (auto batch_sz = 0UL; curr_idx < last_idx && batch_sz < batch_cap;
         ++curr_idx) {
      batch_sz += seqs[curr_idx]->inflated_len;
    }

    return curr_idx;
  };

  // reserve memory for groups
  {
    auto seqs_data_sz = std::transform_reduce(
        seqs.cbegin(), seqs.cend(), 0UL, std::plus<std::size_t>(),
        [](std::unique_ptr<biosoup::NucleicAcid> const& seq) -> std::size_t {
          return seq->inflated_len;
        });

    dst.reserve(std::ceil(static_cast<double>(seqs_data_sz) / target_group_sz));
  }

  // create groups
  {
    for (auto first_idx = 0U; first_idx < seqs.size();) {
      auto const last_idx =
          find_batch_last(first_idx, seqs.size(), target_group_sz);
      dst.push_back(
          OverlapGroup{.first_seq_id = first_idx, .last_seq_id = last_idx});

      first_idx = last_idx;
    }

    for (auto group_id = 0U; group_id < dst.size(); ++group_id) {
      auto const& it = dst[group_id];
      for (auto idx = it.first_seq_id; idx != it.last_seq_id; ++idx) {
        seq_ovlp_group_ids[idx] = group_id;
      }
    }
  }

  auto store_overlap =
      [&dst, &seq_ovlp_group_ids](biosoup::Overlap const& ovlp) -> void {
    dst[seq_ovlp_group_ids[ovlp.lhs_id]].ovlps.push_back(ovlp);
  };

  auto minimizer_engine =
      ram::MinimizerEngine(thread_pool, map_cfg.kmer_len, map_cfg.win_len);

  auto map_futures = std::vector<std::future<std::vector<biosoup::Overlap>>>();

  auto const map_sequence =
      [&minimizer_engine](std::unique_ptr<biosoup::NucleicAcid> const& seq)
      -> std::vector<biosoup::Overlap> {
    return minimizer_engine.Map(seq, true, true, true);
  };

  auto const async_map_sequence =
      [&thread_pool, &map_sequence](
          std::reference_wrapper<std::unique_ptr<biosoup::NucleicAcid> const>
              seq) { return thread_pool->Submit(map_sequence, seq); };

  auto timer = biosoup::Timer();
  timer.Start();

  for (auto const& it : dst) {
    minimizer_engine.Minimize(std::next(seqs.begin(), it.first_seq_id),
                              std::next(seqs.begin(), it.last_seq_id), true);
    minimizer_engine.Filter(map_cfg.filter_p);

    fmt::print(stderr,
               "[camel::FindOverlaps]({:12.3f}) minimized {} / {} sequences\n",
               timer.Stop(), it.last_seq_id, seqs.size());
    timer.Start();

    for (auto map_batch_first = 0U; map_batch_first < it.last_seq_id;) {
      auto const map_batch_last = find_batch_last(
          map_batch_first, it.last_seq_id, detail::kMapBatchCap);

      map_futures.reserve(map_batch_last - map_batch_first);
      std::transform(std::next(seqs.begin(), map_batch_first),
                     std::next(seqs.begin(), map_batch_last),
                     std::back_inserter(map_futures), async_map_sequence);
      for (auto& it : map_futures) {
        for (auto&& ovlp : it.get()) {
          store_overlap(ovlp);
        }
      }

      fmt::print(stderr,
                 "[camel::FindOverlaps]({:12.3f}) mapped {} / {}sequences\n ",
                 timer.Stop(), map_batch_last, it.last_seq_id);
      timer.Start();

      map_batch_first = map_batch_last;
      map_futures.clear();
    }
  }

  auto const n_ovlps = std::transform_reduce(
      dst.cbegin(), dst.cend(), 0UL, std::plus<std::size_t>(),
      [](OverlapGroup const& ovlp_group) -> std::size_t {
        return ovlp_group.ovlps.size();
      });

  fmt::print(stderr, "[camel::FindOverlaps] found {} overlaps\n", n_ovlps);
  return dst;
}

}  // namespace camel
