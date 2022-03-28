#include "camel/coverage.h"

#include <fstream>
#include <iterator>

#include "biosoup/timer.hpp"
#include "cereal/access.hpp"
#include "cereal/archives/binary.hpp"
#include "cereal/specialize.hpp"
#include "cereal/types/string.hpp"
#include "cereal/types/vector.hpp"
#include "edlib.h"
#include "fmt/core.h"
#include "fmt/format.h"

namespace camel {

namespace detail {

static constexpr std::size_t kAlignBatchCap = 1UL << 29UL;  // ~0.5gb

static constexpr std::size_t kDefaultPileStorageFileSz = 1UL << 22UL;  // ~4mb

}  // namespace detail

auto CalculateCoverage(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool, MapCfg const map_cfg,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& seqs)
    -> std::vector<Pile> {
  auto dst = std::vector<Pile>();
  auto ovlps = FindOverlaps(thread_pool, map_cfg, seqs);

  dst.reserve(seqs.size());
  std::transform(seqs.cbegin(), seqs.cend(), std::back_inserter(dst),
                 [](std::unique_ptr<biosoup::NucleicAcid> const& seq) -> Pile {
                   return Pile{.id = seq->id,
                               .seq_name = seq->name,
                               .covgs = std::vector<Coverage>(seq->inflated_len,
                                                              {0, 0, 0, 0})};
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
          ++dst[ovlp.lhs_id].covgs[query_pos].mat;

          ++query_pos;
          ++target_pos;

          break;
        }

        case 1: {  // insertion on target
          ++dst[ovlp.lhs_id].covgs[query_pos].del;
          ++query_pos;
          break;
        }

        case 2: {  // insertion on query
          ++dst[ovlp.lhs_id].covgs[query_pos].ins;
          ++target_id;
          break;
        }

        case 3: {  // mismatch
          ++dst[ovlp.lhs_id].covgs[query_pos].mis;

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

template <class Archive>
auto save(Archive& archive, Coverage const& cov) -> void {
  archive(cov.mat, cov.del, cov.ins, cov.mis);
};

template <class Archive>
auto load(Archive& archive, Coverage& cov) -> void {
  archive(cov.mat, cov.del, cov.ins, cov.mis);
}

template <class Archive>
auto save(Archive& archive, Pile const& pile) -> void {
  archive(pile.id, pile.seq_name, pile.covgs);
}

template <class Archive>
auto load(Archive& archive, Pile& pile) -> void {
  archive(pile.id, pile.seq_name, pile.covgs);
}

auto SerializePiles(std::shared_ptr<thread_pool::ThreadPool> thread_pool,
                    std::vector<Pile> const& piles,
                    std::filesystem::path const& dst_dir) -> void {
  SerializePiles(thread_pool, piles, dst_dir,
                 detail::kDefaultPileStorageFileSz);
}

auto SerializePiles(std::shared_ptr<thread_pool::ThreadPool> thread_pool,
                    std::vector<Pile> const& piles,
                    std::filesystem::path const& dst_dir,
                    std::size_t const expected_file_sz) -> void {
  using PileConstIter = std::vector<Pile>::const_iterator;

  if (std::filesystem::exists(dst_dir)) {
    std::filesystem::remove_all(dst_dir);
  }

  std::filesystem::create_directory(dst_dir);

  auto const find_batch_end = [](PileConstIter const first,
                                 PileConstIter const last,
                                 std::size_t const batch_cap) -> PileConstIter {
    auto curr_iter = first;
    for (auto curr_batch_sz = 0UL;
         curr_batch_sz < batch_cap && curr_iter < last; ++curr_iter) {
      curr_batch_sz += (curr_iter->covgs.size()) * sizeof(Coverage);
    }

    return curr_iter;
  };

  auto const serialize_batch = [&dst_dir](PileConstIter const first,
                                          PileConstIter const last,
                                          std::size_t file_id) -> void {
    auto const dst_path =
        dst_dir / fmt::format("pile_dump_{:04d}.camel", file_id);
    auto dst_fstrm = std::fstream(
        dst_path, std::ios::out | std::ios::trunc | std::ios::binary);

    auto archive = cereal::BinaryOutputArchive(dst_fstrm);

    archive(static_cast<std::size_t>(std::distance(first, last)));
    for (auto it = first; it != last; ++it) {
      archive(*it);
    }
  };

  {
    auto batch_nxt_id = 0U;
    auto ser_futures = std::vector<std::future<void>>();

    for (auto batch_begin = piles.cbegin(); batch_begin != piles.cend();) {
      auto const batch_end =
          find_batch_end(batch_begin, piles.cend(), expected_file_sz);

      ser_futures.emplace_back(thread_pool->Submit(serialize_batch, batch_begin,
                                                   batch_end, batch_nxt_id++));

      batch_begin = batch_end;
    }

    for (auto& it : ser_futures) {
      it.get();
    }
  }
}

auto DeserializePiles(std::shared_ptr<thread_pool::ThreadPool> thread_pool,
                      std::filesystem::path const& src_dir)
    -> std::vector<Pile> {
  auto dst = std::vector<Pile>();

  auto pile_futures = std::vector<std::future<std::vector<Pile>>>();
  for (auto const& it : std::filesystem::directory_iterator(src_dir)) {
    pile_futures.emplace_back(thread_pool->Submit(
        [](std::filesystem::directory_entry const& dir_entry)
            -> std::vector<Pile> {
          auto dst = std::vector<Pile>();

          auto src_fstrm =
              std::fstream(dir_entry.path(), std::ios::in | std::ios::binary);
          auto archive = cereal::BinaryInputArchive(src_fstrm);

          auto sz = std::size_t();
          archive(sz);

          dst.resize(sz);
          for (auto i = 0UL; i < sz; ++i) {
            archive(dst[i]);
          }

          dst.shrink_to_fit();
          return dst;
        },
        it));
  }

  for (auto& it : pile_futures) {
    auto piles = it.get();
    std::move(piles.begin(), piles.end(), std::back_inserter(dst));
  }

  dst.shrink_to_fit();
  return dst;
}

}  // namespace camel
