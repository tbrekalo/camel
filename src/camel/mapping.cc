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
static constexpr auto kOvlpBlockCap = 1UL << 16UL;

}  // namespace detail

MapCfg::MapCfg(std::uint8_t kmer_len, std::uint8_t win_len, double filter_p)
    : kmer_len(kmer_len), win_len(win_len), filter_p(filter_p) {}

auto FindOverlaps(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool, MapCfg const map_cfg,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads)
    -> std::vector<Group> {
  auto dst = std::vector<Group>();

  auto group_ids = std::vector<std::uint32_t>(reads.size());
  auto group_sizes = std::vector<std::uint32_t>(reads.size(), 1U);
  std::iota(group_ids.begin(), group_ids.end(), 0U);

  auto n_ovlps = 0U;
  auto read_ovlp_cnt = std::vector<std::size_t>(reads.size(), 0U);
  auto read_ovlp_vecs = std::vector<std::vector<std::vector<biosoup::Overlap>>>(
      reads.size(), std::vector<std::vector<biosoup::Overlap>>());
  for (auto& it : read_ovlp_vecs) {
    it.emplace_back();
  }

  // auto group_ovlp_cnt = std::vector<std::uint64_t>(reads.size());

  auto minimizer_engine =
      ram::MinimizerEngine(thread_pool, map_cfg.kmer_len, map_cfg.win_len);

  auto map_futures = std::vector<std::future<std::vector<biosoup::Overlap>>>();

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

  auto join_groups = [&group_ids, &group_sizes, &find_group_id](
                         std::uint32_t const lhs_id,
                         std::uint32_t const rhs_id) -> std::uint32_t {
    auto const lhs_group_id = find_group_id(lhs_id);
    auto const rhs_group_id = find_group_id(rhs_id);

    if (lhs_group_id != rhs_group_id) {
      auto& lhs_group_size = group_sizes[lhs_group_id];
      auto& rhs_group_size = group_sizes[rhs_group_id];

      if (lhs_group_size >= rhs_group_size) {
        lhs_group_size += rhs_group_size;
        rhs_group_size = 0U;

        group_ids[rhs_id] = lhs_group_id;
      } else {
        rhs_group_size += lhs_group_size;
        lhs_group_size = 0U;

        group_ids[lhs_id] = rhs_group_id;

        return rhs_group_id;
      }
    }

    return lhs_group_id;
  };

  auto const store_overlap =
      [&read_ovlp_cnt, &read_ovlp_vecs](biosoup::Overlap const& ovlp) -> void {
    auto const query_id = ovlp.lhs_id;
    auto const target_id = ovlp.rhs_id;
    auto& ovlp_vec = read_ovlp_vecs[query_id];
    if (ovlp_vec.back().size() == detail::kOvlpBlockCap) {
      if (ovlp_vec.size() == ovlp_vec.capacity()) {
        ovlp_vec.reserve(ovlp_vec.size() * 1.2);
      }

      ovlp_vec.emplace_back();
    }

    ovlp_vec.back().push_back(ovlp);
    ++read_ovlp_cnt[query_id];
    ++read_ovlp_cnt[target_id];
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
          store_overlap(ovlp);
          ++n_ovlps;
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
  {
    timer.Start();
    auto const n_groups =
        std::count_if(group_sizes.cbegin(), group_sizes.cend(),
                      [](std::uint32_t const sz) -> bool { return sz > 1U; });

    auto group_id_to_handle = tsl::robin_map<std::uint32_t, std::uint32_t>();
    group_id_to_handle.reserve(n_groups);

    for (auto g_id = 0U, g_handle = 0U; g_id < reads.size(); ++g_id) {
      if (group_sizes[g_id] > 1U) {
        group_id_to_handle.emplace(g_id, g_handle++);
      }
    }

    fmt::print(stderr,
               "[camel::FindOverlaps]({:12.3f}) created {} group handles\n",
               timer.Stop(), group_id_to_handle.size());

    timer.Start();
    dst.resize(n_groups);
    for (auto read_id = 0U; read_id < reads.size(); ++read_id) {
      auto g_handle = group_id_to_handle[find_group_id(read_id)];
      auto read_ovlp_vec = std::move(read_ovlp_vecs[read_id]);

      read_ovlp_vec.back().shrink_to_fit();
      std::move(read_ovlp_vec.begin(), read_ovlp_vec.end(),
                std::back_inserter(dst[g_handle].ovlp_vecs));

      if (read_id % 20000U == 0U) {
        fmt::print(stderr, "distrubuted {:12d} / {:12d} read overlap vecs\r",
                   read_id, reads.size());
      }
    }

    auto group_seq_sizes = std::vector<std::uint64_t>(n_groups);
    for (auto read_id = 0U; read_id < reads.size(); ++read_id) {
      auto const g_id = find_group_id(read_id);
      auto const g_handle = group_id_to_handle[g_id];
      group_seq_sizes[g_handle] += reads[read_id]->inflated_len;
    }

    // TODO: plot group information

    for (auto const& it : group_id_to_handle) {
      dst[it.second].read_n_ovlps.reserve(group_sizes[it.first]);
      dst[it.second].ovlp_vecs.shrink_to_fit();
    }

    for (auto read_id = 0U; read_id < reads.size(); ++read_id) {
      auto const g_handle = group_id_to_handle[find_group_id(read_id)];
      dst[g_handle].read_n_ovlps.push_back(ReadIdOvlpCnt{
          .read_id = read_id, .n_overlaps = read_ovlp_cnt[read_id]});
    }

    // for (auto g_id = 0U; g_id < reads.size(); ++g_id) {
    //   if (group_sizes[g_id] > 1) {
    //     fmt::print(stderr, "[de::group_sz] {} : {}\n", group_id_to_handle[g_id],
    //                group_sizes[g_id]);
    //   }
    // }

    auto const invalid_iter = std::remove_if(
        dst.begin(), dst.end(),
        [](Group const& g) -> bool { return g.read_n_ovlps.size() == 1UL; });

    auto const n_invalid = std::distance(invalid_iter, dst.end());

    dst.erase(invalid_iter, dst.end());

    fmt::print(
        stderr,
        "[camel::FindOverlaps]({:12.3f}) transformed {} overlaps in {} "
        "groups\n"
        "[camel::FindOverlaps] largest group is {} bytes\n",
        timer.Stop(), n_ovlps, n_groups - n_invalid,
        *std::max_element(group_seq_sizes.cbegin(), group_seq_sizes.cend()));
  }

  fmt::print(stderr, "[camel::FindOverlaps]({:12.3f})\n", timer.elapsed_time());
  return dst;
}

}  // namespace camel
