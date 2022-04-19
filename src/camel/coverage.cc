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

static constexpr std::size_t kAlignBatchCap = 1UL << 36UL;      // ~68.7gb
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

  std::filesystem::create_directory(pile_storage_dir);
  auto timer = biosoup::Timer();

  // sort groups ascending by read size
  {
    auto group_with_sizes =
        std::vector<std::pair<std::size_t, std::vector<biosoup::Overlap>>>();

    auto const calc_group_sz =
        [&reads](
            std::vector<biosoup::Overlap>::const_iterator first,
            std::vector<biosoup::Overlap>::const_iterator last) -> std::size_t {
      auto unique_read_ids = tsl::robin_set<std::size_t>();
      for (auto it = first; it != last; ++it) {
        unique_read_ids.insert(it->lhs_id);
        unique_read_ids.insert(it->rhs_id);
      }

      return std::transform_reduce(
          unique_read_ids.cbegin(), unique_read_ids.cend(), 0UL,
          std::plus<std::size_t>(),
          [&reads](std::uint32_t const read_id) -> std::size_t {
            return reads[read_id]->inflated_len;
          });
    };

    timer.Start();

    group_with_sizes.reserve(ovlp_groups.size());
    std::transform(
        std::make_move_iterator(ovlp_groups.begin()),
        std::make_move_iterator(ovlp_groups.end()),
        std::back_inserter(group_with_sizes),
        [&calc_group_sz](std::vector<biosoup::Overlap> ovlps)
            -> std::pair<std::size_t, std::vector<biosoup::Overlap>> {
          return std::pair(calc_group_sz(ovlps.cbegin(), ovlps.cend()),
                           std::move(ovlps));
        });

    ovlp_groups.clear();
    std::sort(
        group_with_sizes.begin(), group_with_sizes.end(),
        [](std::pair<std::size_t, std::vector<biosoup::Overlap>> const& lhs,
           std::pair<std::size_t, std::vector<biosoup::Overlap>> const& rhs)
            -> bool { return lhs.first < rhs.first; });

    std::transform(
        std::make_move_iterator(group_with_sizes.begin()),
        std::make_move_iterator(group_with_sizes.end()),
        std::back_inserter(ovlp_groups),
        [](std::pair<std::size_t, std::vector<biosoup::Overlap>> pso)
            -> std::vector<biosoup::Overlap> { return std::move(pso.second); });

    fmt::print(
        stderr,
        "[camel::CalculateCoverage]({:12.3f}) sorted overlap groups by size\n",
        timer.Stop());
  }

  // auto group_sort_futures = std::vector<std::future<void>>();

  // auto const ovlp_cmp_by_lhs_id = [](biosoup::Overlap const& a,
  //                                    biosoup::Overlap const& b) {
  //   return a.lhs_id < b.rhs_id;
  // };

  // auto const sort_ovlp_range =
  //     [](std::vector<biosoup::Overlap>::const_iterator first,
  //        std::vector<biosoup::Overlap>::const_iterator last,
  //        auto&& ovlp_cmp_fn) -> void {
  //   static_assert(
  //       detail::IsOverlapCmpCallableV<decltype(ovlp_cmp_fn)>,
  //       "[camel::FindOverlaps::sort_ovlp_range] ovlp_cmp_fn must be callable"
  //       "with signature (biosoup::Overlap const&, biosoup::Overlap const&) ->
  //       " "bool");

  //   std::sort(first, last, std::forward<decltype(ovlp_cmp_fn)>(ovlp_cmp_fn));
  // };

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

  auto const refresh_active_piles =
      [&reads, &active_piles, &try_init_pile](
          std::vector<biosoup::Overlap>::const_iterator first,
          std::vector<biosoup::Overlap>::const_iterator last) -> void {
    active_piles.clear();
    for (; first != last; ++first) {
      try_init_pile(reads[first->lhs_id]);
      try_init_pile(reads[first->rhs_id]);
    }
  };

  // auto const overlap_strings = 
  //   [&reads](biosoup::Overlap const& ovlp) -> 
  //     std::pair<std::string, std::string> {

  // 
  // };

  for (auto& ovlp_group : ovlp_groups) {
    refresh_active_piles(ovlp_group.cbegin(), ovlp_group.cend());

    ovlp_group.clear();
    ovlp_group.shrink_to_fit();
  }
}

}  // namespace camel
