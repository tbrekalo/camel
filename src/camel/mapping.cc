#include "camel/mapping.h"

#include <algorithm>
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
static constexpr auto kOvlpBuffBlockCap = static_cast<std::size_t>(5e5);

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

  auto ovlp_buff_blocks = std::vector<std::vector<biosoup::Overlap>>();

  ovlp_buff_blocks.emplace_back();
  ovlp_buff_blocks.back().reserve(detail::kOvlpBuffBlockCap);

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
      [&ovlp_buff_blocks](biosoup::Overlap const& ovlp) -> void {
    if (ovlp_buff_blocks.back().size() == detail::kOvlpBuffBlockCap) {
      ovlp_buff_blocks.emplace_back();
      ovlp_buff_blocks.reserve(detail::kOvlpBuffBlockCap);
    }

    ovlp_buff_blocks.back().push_back(ovlp);
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

    for (auto g_id = 0U, h_id = 0U; g_id < reads.size(); ++g_id) {
      if (group_ovlp_cnt[g_id] > 0) {
        group_id_to_handle.emplace(g_id, h_id++);
      }
    }

    group_overlaps.resize(n_groups);
    auto progress_bar = biosoup::ProgressBar(ovlp_buff_blocks.size(), 10U);
    for (auto block_id = 0U; block_id < ovlp_buff_blocks.size(); ++block_id) {
      auto& ovlp_vec = ovlp_buff_blocks[block_id];
      for (auto const& ovlp : ovlp_vec) {
        auto const g_id = find_group_id(ovlp.lhs_id);
        auto const g_handle = group_id_to_handle[g_id];

        if (group_overlaps[g_handle].size() ==
            group_overlaps[g_handle].capacity()) {
          auto const g_ovlp_cnt = group_ovlp_cnt[g_id];
          if (group_overlaps[g_handle].capacity() * 1.55 < g_ovlp_cnt) {
            group_overlaps[g_handle].reserve(group_overlaps.capacity() * 1.5);
          } else {
            group_overlaps[g_handle].reserve(g_ovlp_cnt);
          }
        }

        group_overlaps[g_handle].push_back(ovlp);
      }

      ovlp_vec.clear();
      ovlp_vec.shrink_to_fit();

      ++progress_bar;
      fmt::print("{}\r", progress_bar);
    }

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
        "[camel::FindOverlaps] larges group is {} bytes\n",
        timer.Stop(), n_ovlps, n_groups,
        *std::max_element(group_seq_sizes.cbegin(), group_seq_sizes.cend()));
  }

  fmt::print(stderr, "[camel::FindOverlaps]({:12.3f})\n", timer.elapsed_time());

  return group_overlaps;
}

}  // namespace camel
