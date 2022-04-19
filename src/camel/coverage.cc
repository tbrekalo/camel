#include "camel/coverage.h"

#include <chrono>
#include <deque>
#include <functional>
#include <numeric>

#include "biosoup/timer.hpp"
#include "camel/io.h"
#include "detail/overlap.h"
#include "edlib.h"
#include "fmt/core.h"
#include "tsl/robin_map.h"
#include "tsl/robin_set.h"

namespace camel {

namespace detail {

static constexpr std::size_t kAlignBatchCap = 1UL << 34UL;      // ~17.18gb
static constexpr std::size_t kDefaultSeqGroupSz = 1UL << 32UL;  // 4.3gb

}  // namespace detail

auto CalculateCoverage(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool, MapCfg const map_cfg,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads,
    std::filesystem::path const& pile_storage_dir) -> void {
  auto ovlp_groups = FindOverlaps(thread_pool, map_cfg, reads);

  if (std::filesystem::exists(pile_storage_dir)) {
    std::filesystem::remove_all(pile_storage_dir);
  }

  auto find_batch_last = [&reads](
                             std::vector<ReadIdOvlpCnt>::const_iterator first,
                             std::vector<ReadIdOvlpCnt>::const_iterator last,
                             std::size_t const batch_cap)
      -> std::vector<ReadIdOvlpCnt>::const_iterator {
    for (auto batch_sz = 0UL; first != last && batch_sz < batch_cap; ++first) {
      batch_sz += reads[first->read_id]->inflated_len * sizeof(Coverage) +
                  first->n_overlaps * sizeof(biosoup::Overlap);
    }

    return first;
  };

  std::filesystem::create_directory(pile_storage_dir);
  auto timer = biosoup::Timer();

  // sort groups ascending by read size
  {
    auto group_with_sizes = std::vector<std::pair<std::size_t, Group>>();

    auto const calc_group_sz =
        [&reads](
            std::vector<ReadIdOvlpCnt>::const_iterator first,
            std::vector<ReadIdOvlpCnt>::const_iterator last) -> std::size_t {
      return std::transform_reduce(
          first, last, 0UL, std::plus<std::size_t>(),
          [&reads](ReadIdOvlpCnt r_id_ovlp_cnt) -> std::size_t {
            return reads[r_id_ovlp_cnt.read_id]->inflated_len;
          });
    };

    timer.Start();

    group_with_sizes.reserve(ovlp_groups.size());
    std::transform(
        std::make_move_iterator(ovlp_groups.begin()),
        std::make_move_iterator(ovlp_groups.end()),
        std::back_inserter(group_with_sizes),
        [&calc_group_sz](Group&& group) -> std::pair<std::size_t, Group> {
          return std::pair(calc_group_sz(group.read_n_ovlps.cbegin(),
                                         group.read_n_ovlps.cend()),
                           std::move(group));
        });

    ovlp_groups.clear();
    std::sort(group_with_sizes.begin(), group_with_sizes.end(),
              [](std::pair<std::size_t, Group> const& lhs,
                 std::pair<std::size_t, Group> const& rhs) -> bool {
                return lhs.first < rhs.first;
              });

    std::transform(std::make_move_iterator(group_with_sizes.begin()),
                   std::make_move_iterator(group_with_sizes.end()),
                   std::back_inserter(ovlp_groups),
                   [](std::pair<std::size_t, Group>&& psg) -> Group {
                     return std::move(psg.second);
                   });

    fmt::print(
        stderr,
        "[camel::CalculateCoverage]({:12.3f}) sorted overlap groups by size\n",
        timer.Stop());
  }

  auto active_piles = tsl::robin_map<std::uint32_t, Pile>();

  auto const try_init_pile =
      [&active_piles](
          std::unique_ptr<biosoup::NucleicAcid> const& read) -> void {
    auto& read_pile = active_piles[read->id];
    if (read_pile.seq_name.empty()) {
      read_pile.id = read->id;
      read_pile.seq_name = read->name;
      read_pile.covgs.resize(read->inflated_len);
    }
  };

  auto parallel_gather =
      [&thread_pool](
          ReadIdOvlpCnt const& read_info,
          Group const& active_group) -> std::vector<biosoup::Overlap> {
    auto dst = std::vector<biosoup::Overlap>();
    dst.resize(read_info.n_overlaps);

    auto gather_futures =
        std::vector<std::future<std::vector<biosoup::Overlap>>>();

    auto empty_idx = std::atomic<std::size_t>();
    for (auto const& ovlp_vec : active_group.ovlp_vecs) {
      gather_futures.emplace_back(thread_pool->Submit(
          [&dst, &empty_idx](
              std::uint32_t const read_id,
              std::vector<biosoup::Overlap> const& ovlp_vec) -> void {
            for (auto const& ovlp : ovlp_vec) {
              if (ovlp.lhs_id == read_id) {
                dst[++empty_idx] = ovlp;
              } else if (ovlp.rhs_id == read_id) {
                dst[++empty_idx] = detail::ReverseOverlap(ovlp);
              }
            }
          },
          read_info.read_id, ovlp_vec));
    }

    return dst;
  };


  for (auto g_id = 0U; g_id < ovlp_groups.size(); ++g_id) {
    auto active_group = std::move(ovlp_groups[g_id]);
    for (auto batch_first = active_group.read_n_ovlps.cbegin();
         batch_first != active_group.read_n_ovlps.end();) {
      auto batch_end =
          find_batch_last(batch_first, active_group.read_n_ovlps.cend(),
                          detail::kAlignBatchCap);
    }
  }
}

}  // namespace camel
