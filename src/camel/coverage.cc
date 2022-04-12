#include "camel/coverage.h"

#include <functional>
#include <numeric>

#include "biosoup/timer.hpp"
#include "camel/io.h"
#include "detail/overlap.h"
#include "edlib.h"
#include "fmt/core.h"

namespace camel {

namespace detail {

static constexpr std::size_t kAlignBatchCap = 1UL << 34UL;  // ~16gb

}  // namespace detail

auto CalculateCoverage(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool, MapCfg const map_cfg,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& seqs,
    std::filesystem::path const& pile_storage_dir) -> void {
  auto ovlps = FindOverlaps(thread_pool, map_cfg, seqs);

  if (std::filesystem::exists(pile_storage_dir)) {
    std::filesystem::remove_all(pile_storage_dir);
  }

  auto const find_batch_end =
      [&seqs, &ovlps](std::size_t const begin_idx, std::size_t const end_idx,
                      std::size_t const batch_cap) -> std::size_t {
    auto curr_idx = begin_idx;
    for (auto batch_sz = 0UL; curr_idx < end_idx && batch_sz < batch_cap;
         ++curr_idx) {
      batch_sz += seqs[curr_idx]->inflated_len * sizeof(Coverage);
      batch_sz += std::transform_reduce(
          ovlps[curr_idx].cbegin(), ovlps[curr_idx].cend(), 0UL,
          std::plus<std::size_t>(), detail::OverlapLength);
    }

    return curr_idx;
  };

  auto const overlap_strings =
      [&seqs](
          biosoup::Overlap const& ovlp) -> std::pair<std::string, std::string> {
    auto query_str = seqs[ovlp.lhs_id]->InflateData(
        ovlp.lhs_begin, ovlp.lhs_end - ovlp.lhs_begin);

    auto target_str = seqs[ovlp.rhs_id]->InflateData(
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
    auto dst = Pile{.id = seq->id,
                    .seq_name = seq->name,
                    .covgs = std::vector<Coverage>(seq->inflated_len)};
    for (auto const& ovlp : ovlps) {
      auto [query_substr, target_substr] = overlap_strings(ovlp);
      auto const edlib_res_align = align_strings(query_substr, target_substr);

      query_substr.clear();
      query_substr.shrink_to_fit();

      target_substr.clear();
      target_substr.shrink_to_fit();

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

  auto const async_calc_coverage_for =
      [&thread_pool, &calc_coverage_for](
          std::reference_wrapper<std::unique_ptr<biosoup::NucleicAcid> const>
              seq,
          std::reference_wrapper<std::vector<biosoup::Overlap>> ovlps)
      -> std::future<Pile> {
    return thread_pool->Submit(calc_coverage_for, seq, ovlps);
  };

  {
    auto batch_piles = std::vector<Pile>();
    auto timer = biosoup::Timer();

    auto pile_futures = std::vector<std::future<Pile>>();
    for (auto covg_batch_begin = 0UL, batch_id = 0UL;
         covg_batch_begin < seqs.size(); ++batch_id) {
      timer.Start();
      auto const covg_batch_end =
          find_batch_end(covg_batch_begin, seqs.size(), detail::kAlignBatchCap);

      auto const batch_n_seqs = covg_batch_end - covg_batch_begin;
      pile_futures.reserve(batch_n_seqs);
      batch_piles.reserve(batch_n_seqs);

      std::transform(std::next(seqs.begin(), covg_batch_begin),
                     std::next(seqs.begin(), covg_batch_end),
                     std::next(ovlps.begin(), covg_batch_begin),
                     std::back_inserter(pile_futures), async_calc_coverage_for);

      std::transform(pile_futures.begin(), pile_futures.end(),
                     std::back_inserter(batch_piles),
                     std::mem_fn(&std::future<Pile>::get));

      fmt::print(stderr,
                 "[camel::CalculateCoverage]({:12.3f}) calculated coverage for "
                 "{} reads\n",
                 timer.Stop(), covg_batch_end - covg_batch_begin);

      SerializePileBatch(batch_piles.cbegin(), batch_piles.cend(),
                         pile_storage_dir,
                         fmt::format("pile_batch_{:04d}", batch_id));

      covg_batch_begin = covg_batch_end;
      pile_futures.clear();
      batch_piles.clear();
    }
  }
}

}  // namespace camel
