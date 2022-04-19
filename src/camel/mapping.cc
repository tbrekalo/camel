#include "camel/mapping.h"

#include <algorithm>
#include <deque>
#include <iterator>
#include <numeric>

#include "biosoup/progress_bar.hpp"
#include "biosoup/timer.hpp"
#include "detail/overlap.h"
#include "fmt/core.h"
#include "fmt/ostream.h"
#include "ram/minimizer_engine.hpp"
#include "tsl/robin_map.h"

namespace camel {

namespace detail {

static constexpr auto kMinimizeBatchCap = 1UL << 32UL;
static constexpr auto kMapBatchCap = 1UL << 30UL;
static constexpr auto kOvlpBlockCap = static_cast<std::size_t>(2e6);

}  // namespace detail

MapCfg::MapCfg(std::uint8_t kmer_len, std::uint8_t win_len, double filter_p)
    : kmer_len(kmer_len), win_len(win_len), filter_p(filter_p) {}

auto FindOverlaps(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool, MapCfg const map_cfg,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads)
    -> std::vector<std::vector<biosoup::Overlap>> {
  auto group_ids = std::vector<std::uint32_t>(reads.size());
  auto group_sizes = std::vector<std::uint32_t>(reads.size());
  std::iota(group_ids.begin(), group_ids.end(), 0U);

  auto group_ovlp_cnt = std::vector<std::uint64_t>(reads.size());

  auto minimizer_engine =
      ram::MinimizerEngine(thread_pool, map_cfg.kmer_len, map_cfg.win_len);

  auto map_futures = std::vector<std::future<std::vector<biosoup::Overlap>>>();

  auto ovlp_blocks = std::vector<std::vector<biosoup::Overlap>>();

  ovlp_blocks.emplace_back();
  ovlp_blocks.back().reserve(detail::kOvlpBlockCap);

  auto find_group_id =
      [&group_ids](std::uint32_t const read_id) -> std::uint32_t {
    auto group_id = read_id;
    while (group_id != group_ids[group_id]) {
      group_id = group_ids[group_id];
    }

    for (auto curr_id = read_id; curr_id != group_id;) {
      auto next_id = group_ids[curr_id];
      group_ids[curr_id] = group_id;

      curr_id = next_id;
    }

    return group_id;
  };

  auto join_groups = [&group_ids, &group_sizes, &group_ovlp_cnt,
                      &find_group_id](
                         std::uint32_t const lhs_id,
                         std::uint32_t const rhs_id) -> std::uint32_t {
    auto const lhs_group_id = find_group_id(lhs_id);
    auto const rhs_group_id = find_group_id(rhs_id);

    if (lhs_group_id != rhs_group_id) {
      auto const lhs_group_size = group_sizes[lhs_group_id];
      auto const rhs_group_size = group_sizes[rhs_group_id];

      if (lhs_group_size >= rhs_group_size) {
        group_ids[rhs_group_id] = lhs_group_id;
        group_sizes[lhs_group_id] += rhs_group_size;

        auto& that_cnt = group_ovlp_cnt[rhs_group_id];
        group_ovlp_cnt[lhs_group_id] += that_cnt;
        that_cnt = 0;

        return lhs_group_id;
      } else {
        group_ids[lhs_group_id] = rhs_group_id;
        group_sizes[rhs_group_id] += lhs_group_size;

        auto& that_cnt = group_ovlp_cnt[lhs_group_id];
        group_ovlp_cnt[rhs_group_id] += that_cnt;
        that_cnt = 0;

        return rhs_group_id;
      }
    }

    return lhs_group_id;
  };

  auto const store_overlap =
      [&ovlp_blocks](biosoup::Overlap const& ovlp) -> void {
    if (ovlp_blocks.back().size() == detail::kOvlpBlockCap) {
      if (ovlp_blocks.size() == ovlp_blocks.capacity()) {
        ovlp_blocks.reserve(ovlp_blocks.capacity() * 1.25);
      }

      ovlp_blocks.emplace_back();
      ovlp_blocks.reserve(detail::kOvlpBlockCap);
    }

    ovlp_blocks.back().push_back(ovlp);
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

  auto map_sequence =
      [&minimizer_engine](std::unique_ptr<biosoup::NucleicAcid> const& read)
      -> std::vector<biosoup::Overlap> {
    return minimizer_engine.Map(read, true, true, true);
  };

  auto map_sequence_async =
      [&thread_pool, &map_sequence](
          std::reference_wrapper<std::unique_ptr<biosoup::NucleicAcid> const>
              read) -> std::future<std::vector<biosoup::Overlap>> {
    return thread_pool->Submit(map_sequence, read);
  };

  auto timer = biosoup::Timer();
  for (auto minimize_first = reads.cbegin(); minimize_first != reads.end();) {
    timer.Start();
    auto const minimize_last = find_batch_last(minimize_first, reads.cend(),
                                               detail::kMinimizeBatchCap);

    minimizer_engine.Minimize(minimize_first, minimize_last, true);
    minimizer_engine.Filter(map_cfg.filter_p);

    fmt::print(stderr,
               "[camel::FindOverlaps]({:12.3f}) minimized {} / {} reads\n",
               timer.Stop(), std::distance(reads.cbegin(), minimize_last),
               reads.size());

    for (auto map_first = reads.cbegin(); map_first < minimize_last;) {
      timer.Start();
      auto const map_last =
          find_batch_last(map_first, minimize_last, detail::kMapBatchCap);

      map_futures.reserve(std::distance(map_first, map_last));
      std::transform(map_first, map_last, std::back_inserter(map_futures),
                     map_sequence_async);

      for (auto& it : map_futures) {
        for (auto const& ovlp : it.get()) {
          auto const g_id = join_groups(ovlp.lhs_id, ovlp.rhs_id);
          ++group_ovlp_cnt[g_id];
          store_overlap(ovlp);
        }
      }

      fmt::print(stderr,
                 "[camel::FindOverlaps]({:12.3f}) mapped {} / {} reads\n",
                 timer.Stop(), std::distance(reads.cbegin(), map_last),
                 std::distance(reads.cbegin(), minimize_last));

      map_first = map_last;
      map_futures.clear();
    }

    minimize_first = minimize_last;
  }

  // group overlaps in groups
  auto group_overlaps = std::vector<std::vector<biosoup::Overlap>>();
  {
    timer.Start();

    auto const [n_groups, n_ovlps] = std::transform_reduce(
        group_ovlp_cnt.cbegin(), group_ovlp_cnt.cend(), std::pair(0U, 0UL),
        [](auto const& lhs, auto const& rhs) {
          return std::pair(lhs.first + rhs.first, lhs.second + rhs.second);
        },
        [](std::uint64_t const cnt) { return std::pair(cnt > 0, cnt); });

    auto group_id_to_handle = tsl::robin_map<std::uint32_t, std::uint32_t>();
    group_id_to_handle.reserve(n_groups);

    for (auto read_id = 0U, h_id = 0U; read_id < reads.size(); ++read_id) {
      auto const g_id = find_group_id(read_id);
      auto const& g_handle = group_id_to_handle[g_id];
      if (g_handle == 0) {
        group_id_to_handle.emplace(g_id, h_id++);
      }
    }

    fmt::print(stderr,
               "[camel::FindOverlaps]({:12.3f}) created {} group handles\n",
               timer.Stop(), group_id_to_handle.size());

    timer.Start();
    auto sort_futures =
        std::vector<std::future<std::vector<std::vector<biosoup::Overlap>>>>();
    sort_futures.reserve(ovlp_blocks.size());

    // initial groupings
    for (auto& ovlp_block : ovlp_blocks) {
      sort_futures.emplace_back(thread_pool->Submit(
          [&group_ids, &group_id_to_handle,
           ovlp_vec = std::move(ovlp_block)](std::size_t const n_groups) mutable
          -> std::vector<std::vector<biosoup::Overlap>> {
            auto dst = std::vector<std::vector<biosoup::Overlap>>(
                n_groups, std::vector<biosoup::Overlap>());
            for (auto const& it : ovlp_vec) {
              dst[group_id_to_handle[group_ids[it.lhs_id]]].push_back(it);
            }

            ovlp_vec.clear();
            ovlp_vec.shrink_to_fit();

            return dst;
          },
          n_groups));
    }

    auto task_pack_queue =
        std::deque<std::vector<std::vector<biosoup::Overlap>>>();

    for (auto& it : sort_futures) {
      task_pack_queue.push_back(it.get());
    }

    fmt::print(
        stderr,
        "[camel::FindOverlaps]({:12.3f}) sorted initial {} overlap blocks\n",
        timer.Stop(), task_pack_queue.size());

    timer.Start();

    sort_futures.clear();
    while (true) {
      while (task_pack_queue.size() >= 2) {
        auto lhs_vec = std::move(task_pack_queue.front());
        task_pack_queue.pop_front();

        auto rhs_vec = std::move(task_pack_queue.front());
        task_pack_queue.pop_front();

        sort_futures.emplace_back(
            thread_pool->Submit([lhs_ovlps = std::move(lhs_vec),
                                 rhs_ovlps = std::move(rhs_vec)]() mutable
                                -> std::vector<std::vector<biosoup::Overlap>> {
              auto dst = std::vector<std::vector<biosoup::Overlap>>();
              dst.reserve(lhs_ovlps.size());

              for (auto idx = 0U; idx < lhs_ovlps.size(); ++idx) {
                if (lhs_ovlps[idx].size() < rhs_ovlps[idx].size()) {
                  std::swap(lhs_ovlps[idx], rhs_ovlps[idx]);
                }

                lhs_ovlps[idx].reserve(lhs_ovlps[idx].size() +
                                       rhs_ovlps[idx].size());

                std::copy(rhs_ovlps[idx].cbegin(), rhs_ovlps[idx].cend(),
                          std::back_inserter(lhs_ovlps[idx]));

                dst.push_back(std::move(lhs_ovlps[idx]));

                rhs_ovlps[idx].clear();
                rhs_ovlps[idx].shrink_to_fit();
              }

              lhs_ovlps.clear();
              lhs_ovlps.shrink_to_fit();

              rhs_ovlps.clear();
              rhs_ovlps.shrink_to_fit();

              return dst;
            }));
      }

      for (auto& it : sort_futures) {
        task_pack_queue.push_back(it.get());
      }

      sort_futures.clear();

      if (task_pack_queue.size() == 1UL) {
        group_overlaps = std::move(task_pack_queue.front());
        break;
      }
    }

    // group_overlaps.resize(n_groups);
    // auto progress_bar = biosoup::ProgressBar(ovlp_buff_blocks.size(), 10U);
    // for (auto block_id = 0U; block_id < ovlp_buff_blocks.size(); ++block_id)
    // {
    //   auto& ovlp_vec = ovlp_buff_blocks[block_id];
    //   for (auto const& ovlp : ovlp_vec) {
    //     auto const g_id = find_group_id(ovlp.lhs_id);
    //     auto const g_handle = group_id_to_handle[g_id];

    //     if (group_overlaps[g_handle].size() ==
    //         group_overlaps[g_handle].capacity()) {
    //       auto const g_ovlp_cnt = group_ovlp_cnt[g_id];
    //       if (group_overlaps[g_handle].capacity() * 1.26 < g_ovlp_cnt) {
    //         group_overlaps[g_handle].reserve(group_overlaps.capacity()
    //         * 1.25);
    //       } else {
    //         group_overlaps[g_handle].reserve(g_ovlp_cnt);
    //       }
    //     }

    //     group_overlaps[g_handle].push_back(ovlp);
    //   }

    //   ovlp_vec.clear();
    //   ovlp_vec.shrink_to_fit();

    //   ++progress_bar;
    //   fmt::print("{}\r", progress_bar);
    // }

    auto group_seq_sizes = std::vector<std::uint64_t>(n_groups);
    for (auto read_id = 0U; read_id < reads.size(); ++read_id) {
      auto const g_id = find_group_id(read_id);
      auto const g_handle = group_id_to_handle[g_id];
      group_seq_sizes[g_handle] += reads[read_id]->inflated_len;
    }

    fmt::print(
        stderr,
        "[camel::FindOverlaps]({:12.3f}) transformed {} overlaps in {} "
        "groups\n"
        "[camel::FindOverlaps] largest group is {} bytes\n",
        timer.Stop(), n_ovlps, n_groups,
        *std::max_element(group_seq_sizes.cbegin(), group_seq_sizes.cend()));
  }

  fmt::print(stderr, "[camel::FindOverlaps]({:12.3f})\n", timer.elapsed_time());
  return group_overlaps;
}

}  // namespace camel
