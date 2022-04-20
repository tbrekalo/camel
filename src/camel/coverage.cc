#include "camel/coverage.h"

#include <chrono>
#include <deque>
#include <functional>
#include <numeric>

#include "biosoup/timer.hpp"
#include "camel/io.h"
#include "detail/overlap.h"
#include "edlib.h"
#include "fmt/compile.h"
#include "fmt/core.h"
#include "tsl/robin_map.h"

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
      batch_sz += reads[first->read_id]->inflated_len * sizeof(Coverage);
      // + first->n_overlaps * sizeof(biosoup::Overlap);
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
          [&reads](ReadIdOvlpCnt ri_oc) -> std::size_t {
            return reads[ri_oc.read_id]->inflated_len;
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

  auto read_id_to_pile_idx = tsl::robin_map<std::uint32_t, std::uint32_t>();
  auto active_piles = std::vector<Pile>();

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

    return std::pair(std::move(query_str), std::move(target_str));
  };

  auto const align_strings =
      [](std::string const& query_str,
         std::string const& target_str) -> EdlibAlignResult {
    return edlibAlign(
        query_str.c_str(), query_str.size(), target_str.c_str(),
        target_str.size(),
        edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));
  };

  auto const update_coverage =
      [&reads, &read_id_to_pile_idx, &active_piles, &overlap_strings,
       &align_strings](biosoup::Overlap const& ovlp) -> void {
    auto [query_substr, target_substr] = overlap_strings(ovlp);
    auto const edlib_res_align = align_strings(query_substr, target_substr);

    query_substr.clear();
    query_substr.shrink_to_fit();

    auto query_id = ovlp.lhs_id;
    auto target_id = ovlp.rhs_id;

    auto query_pos = ovlp.lhs_begin;
    auto target_pos = 0U;

    auto const pile_idx = read_id_to_pile_idx[query_id];
    auto& pile = active_piles[pile_idx];

    for (auto i = 0U; i < edlib_res_align.alignmentLength; ++i) {
      switch (edlib_res_align.alignment[i]) {
        case 0:
        case 3: {  // match
          /* clang-format off */
            switch (target_substr[target_pos]) {
              case 'A': ++pile.covgs[query_pos].a; break;
              case 'C': ++pile.covgs[query_pos].c; break;
              case 'G': ++pile.covgs[query_pos].g; break;
              case 'T': ++pile.covgs[query_pos].t; break;
              default: break;
            }
          /* clang-format on */

          ++query_pos;
          ++target_pos;

          break;
        }

        case 1: {  // insertion on target
          ++pile.covgs[query_pos].del;
          ++query_pos;
          break;
        }

        case 2: {  // insertion on query
          ++pile.covgs[query_pos].ins;
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

  auto const search_and_update =
      [&update_coverage](
          ReadIdOvlpCnt ri_oc,
          std::vector<std::vector<biosoup::Overlap>>::const_iterator first,
          std::vector<std::vector<biosoup::Overlap>>::const_iterator last)
      -> void {
    auto cnt = 0UL;
    for (; first != last; ++first) {
      for (auto const& ovlp : *first) {
        if (ovlp.lhs_id == ri_oc.read_id) {
          update_coverage(ovlp);
          if (++cnt == ri_oc.n_overlaps) {
            return;
          }
        } else if (ovlp.rhs_id == ri_oc.read_id) {
          update_coverage(detail::ReverseOverlap(ovlp));
          if (++cnt == ri_oc.n_overlaps) {
            break;
          }
        }
      }
    }
  };

  auto coverage_futures = std::vector<std::future<void>>();
  auto serialize_futures = std::deque<std::future<
      std::tuple<std::filesystem::path, std::size_t, std::size_t, double>>>();

  auto is_ser_future_ready =
      [](decltype(serialize_futures)::const_reference f) -> bool {
    return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
  };

  auto pop_and_report = [&serialize_futures]() -> void {
    auto const [dst_file, raw_sz, com_sz, comp_time] =
        serialize_futures.front().get();
    serialize_futures.pop_front();

    fmt::print(stderr,
               FMT_COMPILE("[camel::CalculateCoverage]({:12.3f}) {} -> "
                           "comp ratio: {:1.3f}\n"),
               comp_time, dst_file.string(), 1. * raw_sz / com_sz);
  };

  for (auto g_id = 0U; g_id < ovlp_groups.size(); ++g_id) {
    auto active_group = std::move(ovlp_groups[g_id]);
    auto batch_id = 0U;
    for (auto batch_first = active_group.read_n_ovlps.cbegin();
         batch_first != active_group.read_n_ovlps.end(); ++batch_id) {
      timer.Start();
      auto batch_last =
          find_batch_last(batch_first, active_group.read_n_ovlps.cend(),
                          detail::kAlignBatchCap);

      auto const batch_size = std::distance(batch_first, batch_last);
      read_id_to_pile_idx.reserve(batch_size);
      coverage_futures.reserve(batch_size);
      active_piles.reserve(batch_size);

      while (!serialize_futures.empty() &&
             is_ser_future_ready(serialize_futures.front())) {
        pop_and_report();
      }

      for (auto it = batch_first; it != batch_last; ++it) {
        read_id_to_pile_idx.emplace(it->read_id, active_piles.size());
        active_piles.push_back(Pile{
            .id = it->read_id,
            .seq_name = reads[it->read_id]->name,
            .covgs = std::vector<Coverage>(reads[it->read_id]->inflated_len)});
      }

      fmt::print(stderr,
                 "[camel::CalculateCoverage]({:12.3f}) initialized {} piles\n",
                 timer.Stop(), batch_size);
      timer.Start();

      for (auto it = batch_first; it != batch_last; ++it) {
        if (it->n_overlaps > 0) {
          coverage_futures.emplace_back(thread_pool->Submit(
              search_and_update, *it, active_group.ovlp_vecs.cbegin(),
              active_group.ovlp_vecs.cend()));
        }
      }

      std::for_each(coverage_futures.begin(), coverage_futures.end(),
                    std::mem_fn(&std::future<void>::wait));
      coverage_futures.clear();

      fmt::print(stderr,
                 "[camel::CalculateCoverage]({:12.3f}) calculated coverage for "
                 "{} reads\n",
                 timer.Stop(), batch_size);

      serialize_futures.emplace_back(thread_pool->Submit(
          [&pile_storage_dir, group_id = g_id,
           batch_piles =
               std::move(active_piles)](std::uint32_t const batch_id) mutable
          -> std::tuple<std::filesystem::path, std::size_t, std::size_t,
                        double> {
            auto ser_timer = biosoup::Timer();

            ser_timer.Start();
            auto dst_file = SerializePileBatch(
                batch_piles.cbegin(), batch_piles.cend(), pile_storage_dir,
                fmt::format(FMT_COMPILE("group_{:04d}_pile_batch_{:04d}"),
                            group_id, batch_id));

            auto const uncompressed_bytes = std::transform_reduce(
                std::make_move_iterator(batch_piles.begin()),
                std::make_move_iterator(batch_piles.end()), 0UL,
                std::plus<std::size_t>(), [](Pile&& pile) -> std::size_t {
                  return sizeof(Pile::id) + pile.seq_name.size() +
                         pile.covgs.size() * sizeof(Coverage);
                });

            auto const compressed_bytes = std::filesystem::file_size(dst_file);

            return std::tuple(std::move(dst_file), uncompressed_bytes,
                              compressed_bytes, ser_timer.Stop());

            return {};
          },
          batch_id));

      read_id_to_pile_idx.clear();
      active_piles.clear();
      batch_first = batch_last;
    }
  }

  while (!serialize_futures.empty()) {
    pop_and_report();
  }
}

}  // namespace camel
