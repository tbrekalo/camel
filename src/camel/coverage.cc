#include "camel/coverage.h"

#include "biosoup/timer.hpp"
#include "edlib.h"
#include "fmt/core.h"

namespace camel {

namespace detail {

static constexpr std::size_t kAlignBatchCap = 1UL << 29UL;

}

auto CalculateCoverage(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool, MapCfg const map_cfg,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& seqs)
    -> std::vector<std::vector<Coverage>> {
  auto dst = std::vector<std::vector<Coverage>>();
  auto ovlps = FindOverlaps(thread_pool, map_cfg, seqs);

  dst.reserve(seqs.size());
  std::transform(
      seqs.cbegin(), seqs.cend(), std::back_inserter(dst),
      [](std::unique_ptr<biosoup::NucleicAcid> const& seq)
          -> std::vector<Coverage> {
        return std::vector<Coverage>(seq->inflated_len, {0, 0, 0, 0});
      });

  auto const find_batch_end =
      [&seqs](std::size_t const begin_idx, std::size_t const end_idx,
              std::size_t const batch_cap) -> std::size_t {
    auto curr_idx = begin_idx;
    for (auto batch_sz = 0UL; curr_idx < end_idx && batch_sz < batch_cap;
         ++curr_idx) {
      batch_sz += seqs[curr_idx]->inflated_len;
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

  auto const update_coverage = [&overlap_strings, &align_strings,
                                &dst](biosoup::Overlap const& ovlp) -> void {
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
        case 0: {  // match
          ++dst[ovlp.lhs_id][query_pos].mat;

          ++query_pos;
          ++target_pos;

          break;
        }

        case 1: {  // insertion on target
          ++dst[ovlp.lhs_id][query_pos].del;
          ++query_pos;
          break;
        }

        case 2: {  // insertion on query
          ++dst[ovlp.lhs_id][query_pos].ins;
          ++target_id;
          break;
        }

        case 3: {  // mismatch
          ++dst[ovlp.lhs_id][query_pos].mis;

          ++query_pos;
          ++target_pos;
          break;
        }
        default: {
          break;
        }
      }
    }

    edlibFreeAlignResult(edlib_res_align);
  };

  auto const calc_coverage_for =
      [&update_coverage](std::vector<biosoup::Overlap> const& ovlps) -> void {
    std::for_each(ovlps.cbegin(), ovlps.cend(), update_coverage);
  };

  auto const async_calc_coverage_for =
      [&thread_pool, &calc_coverage_for](
          std::reference_wrapper<std::vector<biosoup::Overlap> const> ovlps)
      -> std::future<void> {
    return thread_pool->Submit(calc_coverage_for, ovlps);
  };

  {
    auto timer = biosoup::Timer();

    auto align_futures = std::vector<std::future<void>>();
    for (auto align_batch_begin = 0UL; align_batch_begin < seqs.size();) {
      timer.Start();
      auto const align_batch_end = find_batch_end(
          align_batch_begin, seqs.size(), detail::kAlignBatchCap);

      align_futures.reserve(align_batch_end - align_batch_begin);

      std::transform(std::next(ovlps.begin(), align_batch_begin),
                     std::next(ovlps.begin(), align_batch_end),
                     std::back_inserter(align_futures),
                     async_calc_coverage_for);

      for (auto& it : align_futures) {
        it.wait();
      }

      fmt::print(stderr,
                 "[camel::CalculateCoverage]({:12.3f}) calculated coverage for "
                 "{} reads\n",
                 timer.Stop(), align_batch_end - align_batch_begin);

      align_batch_begin = align_batch_end;
      align_futures.clear();
    }
  }

  return dst;
}

}  // namespace camel
