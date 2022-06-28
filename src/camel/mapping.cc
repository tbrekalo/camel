#include "camel/mapping.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <deque>
#include <iterator>
#include <numeric>

#include "biosoup/timer.hpp"
#include "detail/overlap.h"
#include "fmt/core.h"
#include "ram/minimizer_engine.hpp"
#include "tsl/robin_set.h"

namespace camel {

namespace detail {

static constexpr auto kMinimizeBatchCap = 1UL << 34UL;
static constexpr auto kMapBatchCap = 1UL << 26UL;
static constexpr auto kMxOvlps = 4U * 16U;

static auto FindBatchLast(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator first,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator last,
    std::uint64_t batch_cap)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator {
  for (auto batch_sz = 0UL; batch_sz < batch_cap && first != last; ++first) {
    batch_sz += (*first)->inflated_len;
  }

  return first;
}

}  // namespace detail

MapCfg::MapCfg(std::uint8_t kmer_len, std::uint8_t win_len, double filter_p)
    : kmer_len(kmer_len), win_len(win_len), filter_p(filter_p) {}

auto FindOverlaps(
    State& state, MapCfg const map_cfg,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads)
    -> std::vector<std::vector<biosoup::Overlap>> {
  auto dst = std::vector<std::vector<biosoup::Overlap>>(reads.size());
  auto is_contained = std::vector<std::uint8_t>(reads.size(), 0U);

  auto const mark_containment =
      [&dst, &reads, &is_contained](biosoup::Overlap const& ovlp) -> void {
    auto const kOvlpType =
        detail::DetermineOverlapType(ovlp, reads[ovlp.lhs_id]->inflated_len,
                                     reads[ovlp.rhs_id]->inflated_len);

    switch (kOvlpType) {
      case detail::OverlapType::kLhsContained: {
        is_contained[ovlp.lhs_id] = 1U;
        std::vector<biosoup::Overlap>().swap(dst[ovlp.lhs_id]);
        break;
      }
      case detail::OverlapType::kRhsContained: {
        is_contained[ovlp.rhs_id] = 1U;
        std::vector<biosoup::Overlap>().swap(dst[ovlp.rhs_id]);
      }
      case detail::OverlapType::kInternal: {
        if (1.1 * reads[ovlp.lhs_id]->inflated_len <
            reads[ovlp.rhs_id]->inflated_len) {
          is_contained[ovlp.lhs_id] = 1U;
          std::vector<biosoup::Overlap>().swap(dst[ovlp.lhs_id]);
        } else if (reads[ovlp.lhs_id]->inflated_len >
                   1.1 * reads[ovlp.lhs_id]->inflated_len) {
          is_contained[ovlp.rhs_id] = 1U;
          std::vector<biosoup::Overlap>().swap(dst[ovlp.rhs_id]);
        }
      }
      default: {
        break;
      }
    }
  };

  auto const store_ovlp =
      [&dst, &reads, &is_contained](biosoup::Overlap const& ovlp) -> void {
    if (!is_contained[ovlp.lhs_id]) {
      dst[ovlp.lhs_id].push_back(ovlp);
    }

    if (!is_contained[ovlp.rhs_id]) {
      dst[ovlp.rhs_id].push_back(detail::ReverseOverlap(ovlp));
    }
  };

  auto minimizer_engine = ram::MinimizerEngine(
      state.thread_pool, map_cfg.kmer_len, map_cfg.win_len);

  auto const ovlp_len = [](biosoup::Overlap const& ovlp) -> std::uint32_t {
    return std::max(ovlp.lhs_end - ovlp.lhs_begin,
                    ovlp.rhs_end - ovlp.lhs_begin);
  };

  auto const ovlp_err = [&ovlp_len](biosoup::Overlap const& ovlp) -> double {
    return 1.0 - static_cast<double>(std::min(ovlp.lhs_end - ovlp.lhs_begin,
                                              ovlp.rhs_end - ovlp.rhs_begin)) /
                     static_cast<double>(ovlp_len(ovlp));
  };

  auto const find_batch_last =
      [](std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator
             first,
         std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator
             last,
         std::uint64_t const kBatchCap)
      -> std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator {
    for (auto batch_sz = 0UL; batch_sz < kBatchCap && first != last; ++first) {
      batch_sz += (*first)->inflated_len;
    }

    return first;
  };

  auto timer = biosoup::Timer();

  {
    auto map_futures =
        std::vector<std::future<std::vector<biosoup::Overlap>>>();
    for (auto minimize_first = reads.cbegin();
         minimize_first != reads.cend();) {
      timer.Start();
      auto const minimize_last = find_batch_last(minimize_first, reads.cend(),
                                                 detail::kMinimizeBatchCap);

      minimizer_engine.Minimize(minimize_first, minimize_last, true);
      fmt::print(stderr,
                 "[camel::FindOverlaps]({:12.3f}) minimized "
                 "{} / {} reads\n",
                 timer.Stop(), std::distance(reads.cbegin(), minimize_last),
                 reads.size());

      timer.Start();
      for (auto map_first = minimize_first; map_first != minimize_last;) {
        auto const map_last =
            find_batch_last(map_first, minimize_last, detail::kMapBatchCap);

        map_futures.reserve(std::distance(map_first, map_last));
        for (auto iter = map_first; iter != map_last; ++iter) {
          map_futures.emplace_back(state.thread_pool->Submit(
              [&reads, &minimizer_engine, &ovlp_len,
               &ovlp_err](std::vector<
                          std::unique_ptr<biosoup::NucleicAcid>>::const_iterator
                              read_iter) -> std::vector<biosoup::Overlap> {
                auto ovlps = minimizer_engine.Map(*read_iter, true, true);
                ovlps.erase(
                    std::remove_if(
                        ovlps.begin(), ovlps.end(),
                        [&ovlp_err](biosoup::Overlap const& ovlp) -> bool {
                          return ovlp_err(ovlp) > 0.2;
                        }),

                    ovlps.end());

                if (ovlps.size() > detail::kMxOvlps) {
                  std::nth_element(
                      ovlps.begin(), ovlps.begin() + detail::kMxOvlps,
                      ovlps.end(),
                      [&ovlp_len](biosoup::Overlap const& a,
                                  biosoup::Overlap const& b) -> bool {
                        return ovlp_len(a) > ovlp_len(b);
                      });
                  decltype(ovlps)(ovlps.begin(),
                                  ovlps.begin() + detail::kMxOvlps)
                      .swap(ovlps);

                  return ovlps;
                }

                return ovlps;
              },
              iter));
        }

        for (auto& future : map_futures) {
          auto&& ovlps = future.get();
          for (auto const& ovlp : ovlps) {
            mark_containment(ovlp);
            store_ovlp(ovlp);
          }
        }

        fmt::print(stderr,
                   "\r[camel::FindOverlaps]({:12.3f}) mapped {} / {} reads",
                   timer.Lap(), std::distance(minimize_first, map_last),
                   std::distance(minimize_first, minimize_last));

        map_futures.clear();
        map_first = map_last;
      }

      timer.Stop();
      fmt::print(stderr, "\n");

      minimize_first = minimize_last;
    }
  }

  {
    timer.Start();
    auto defrag_futures = std::vector<std::future<void>>();

    for (auto read_id = 0U; read_id < reads.size(); ++read_id) {
      if (dst[read_id].size() > detail::kMxOvlps) {
        defrag_futures.emplace_back(state.thread_pool->Submit(
            [&dst, &ovlp_len](std::uint32_t const read_id) -> void {
              std::nth_element(dst[read_id].begin(),
                               dst[read_id].begin() + detail::kMxOvlps,
                               dst[read_id].end(),
                               [&ovlp_len](biosoup::Overlap const& a,
                                           biosoup::Overlap const& b) -> bool {
                                 return ovlp_len(a) > ovlp_len(b);
                               });

              std::vector<biosoup::Overlap>(
                  dst[read_id].begin(), dst[read_id].begin() + detail::kMxOvlps)
                  .swap(dst[read_id]);
            },
            read_id));
      }
    }

    std::for_each(defrag_futures.begin(), defrag_futures.end(),
                  std::mem_fn(&std::future<void>::wait));

    fmt::print(
        stderr,
        "[camel::FindOverlaps]({:12.3f}) discarded low quality overlaps\n",
        timer.Stop());
  }

  auto const kNContained =
      std::count_if(dst.begin(), dst.end(),
                    std::mem_fn(&std::vector<biosoup::Overlap>::empty));
  auto const kNTargets = dst.size() - kNContained;

  fmt::print(stderr,
             "[camel::FindOverlaps]({:12.3f}) (kNTargets, kNContained) :="
             " ({}, {})\n",
             timer.elapsed_time(), kNTargets, kNContained);

  return dst;
}

}  // namespace camel
