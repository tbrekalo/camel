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

  auto parallel_gather =
      [&thread_pool](
          ReadIdOvlpCnt const& read_info,
          Group const& active_group) -> std::vector<biosoup::Overlap> {
    auto dst = std::vector<biosoup::Overlap>();
    dst.resize(read_info.n_overlaps);

    auto gather_futures = std::vector<std::future<void>>();

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

  auto const overlap_strings =
      [&reads](
          biosoup::Overlap const& ovlp) -> std::pair<std::string, std::string> {
    auto query_str = reads[ovlp.lhs_id]->InflateData(
        ovlp.lhs_begin, ovlp.lhs_end - ovlp.lhs_begin);

    auto target_str = reads[ovlp.rhs_id]->InflateData(
        ovlp.rhs_begin, ovlp.rhs_end - ovlp.rhs_begin);

    if (!ovlp.strand) {
      auto rc = biosoup::NucleicAcid("", target_str);
      rc.ReverseAndComplement();

      target_str = rc.InflateData();
    }

    return std::make_pair(std::move(query_str), std::move(target_str));
  };

  auto const align_strings =
      [](std::string const& query_str,
         std::string const& target_str) -> EdlibAlignResult {
    return edlibAlign(
        query_str.c_str(), query_str.size(), target_str.c_str(),
        target_str.size(),
        edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));
  };

  auto const calc_coverage_for =
      [&overlap_strings, &align_strings](
          std::unique_ptr<biosoup::NucleicAcid> const& seq,
          std::vector<biosoup::Overlap>& ovlps) -> Pile {
    auto dst =
        Pile{.id = seq->id,
             .seq_name = seq->name,
             .covgs = std::vector<Coverage>(
                 seq->inflated_len,
                 Coverage{.a = 0, .c = 0, .g = 0, .t = 0, .del = 0, .ins = 0})};

    for (auto const& ovlp : ovlps) {
      auto [query_substr, target_substr] = overlap_strings(ovlp);
      auto const edlib_res_align = align_strings(query_substr, target_substr);

      query_substr.clear();
      query_substr.shrink_to_fit();

      auto query_id = ovlp.lhs_id;
      auto target_id = ovlp.rhs_id;

      auto query_pos = ovlp.lhs_begin;
      auto target_pos = 0U;

      for (auto i = 0U; i < edlib_res_align.alignmentLength; ++i) {
        switch (edlib_res_align.alignment[i]) {
          case 0:
          case 3: {  // match
            /* clang-format off */
            switch (target_substr[target_pos]) {
              case 'A': ++dst.covgs[query_pos].a; break;
              case 'C': ++dst.covgs[query_pos].c; break;
              case 'G': ++dst.covgs[query_pos].g; break;
              case 'T': ++dst.covgs[query_pos].t; break;
              default: break;
            }
            /* clang-format on */

            ++query_pos;
            ++target_pos;

            break;
          }

          case 1: {  // insertion on target
            ++dst.covgs[query_pos].del;
            ++query_pos;
            break;
          }

          case 2: {  // insertion on query
            ++dst.covgs[query_pos].ins;
            ++target_id;
            break;
          }

          default: {
            break;
          }
        }
      }

      edlibFreeAlignResult(edlib_res_align);
    };

    ovlps.clear();
    ovlps.shrink_to_fit();

    return dst;
  };


  for (auto g_id = 0U; g_id < ovlp_groups.size(); ++g_id) {
    auto active_group = std::move(ovlp_groups[g_id]);

    for (auto batch_first = active_group.read_n_ovlps.cbegin();
         batch_first != active_group.read_n_ovlps.end();) {
      auto batch_end =
          find_batch_last(batch_first, active_group.read_n_ovlps.cend(),
                          detail::kAlignBatchCap);

      auto overlaps = parallel_gather(*batch_first, active_group);
      auto piles = calc_coverage_for(reads[batch_first->read_id], overlaps);

      active_piles.insert({batch_first->read_id, std::move(piles)});
    }
  }
}

}  // namespace camel
