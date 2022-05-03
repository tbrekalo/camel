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
#include "fmt/ostream.h"
#include "tsl/robin_map.h"

namespace camel {

namespace detail {

static constexpr std::size_t kAlignBatchCap = 1UL << 32UL;      // ~4gb
static constexpr std::size_t kDefaultSeqGroupSz = 1UL << 32UL;  // ~4.3gb

struct InsSignal {
  using ValueType = std::uint16_t;

  ValueType a;
  ValueType c;
  ValueType t;
  ValueType g;
};

[[nodiscard]] static auto CmpReadsOvlpsPairsByReadId(
    ReadOverlapsPair const& lhs, ReadOverlapsPair const& rhs) -> bool {
  return lhs.read->id < rhs.read->id;
}

[[nodiscard]] static auto FindBatchLast(
    std::vector<ReadOverlapsPair>::const_iterator first,
    std::vector<ReadOverlapsPair>::const_iterator last, std::size_t batch_cap) {
  for (auto batch_sz = 0UL; first != last && batch_sz < batch_cap; ++first) {
    batch_sz += first->read->inflated_len * sizeof(Coverage);
  }

  return first;
}

[[nodiscard]] static auto OverlapStrings(
    std::vector<ReadOverlapsPair> const& reads_overlaps,
    biosoup::Overlap const& ovlp) -> std::pair<std::string, std::string> {
  auto const query_read_iter = std::lower_bound(
      reads_overlaps.cbegin(), reads_overlaps.cend(), ovlp.lhs_id,
      [](ReadOverlapsPair const& ro, std::uint32_t const query_id) -> bool {
        return ro.read->id < query_id;
      });

  auto const target_read_iter = std::lower_bound(
      reads_overlaps.cbegin(), reads_overlaps.cend(), ovlp.rhs_id,
      [](ReadOverlapsPair const& ro, std::uint32_t const target_id) -> bool {
        return ro.read->id < target_id;
      });

  auto query_str = query_read_iter->read->InflateData(
      ovlp.lhs_begin, ovlp.lhs_end - ovlp.lhs_begin);

  auto target_str = target_read_iter->read->InflateData(
      ovlp.rhs_begin, ovlp.rhs_end - ovlp.rhs_begin);

  if (!ovlp.strand) {
    auto rc = biosoup::NucleicAcid("", target_str);
    rc.ReverseAndComplement();

    target_str = rc.InflateData();
  }

  return std::pair(std::move(query_str), std::move(target_str));
}

[[nodiscard]] static auto AlignStrings(std::string const& query_str,
                                       std::string const& target_str)
    -> EdlibAlignResult {
  return edlibAlign(
      query_str.c_str(), query_str.size(), target_str.c_str(),
      target_str.size(),
      edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));
}

}  // namespace detail

CAMEL_EXPORT auto CalculateCoverage(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::vector<ReadOverlapsPair> const& reads_overlaps,
    std::filesystem::path const& pile_storage_dir) -> void {
  // consider making reads_overlaps just &; not const&
  if (!std::is_sorted(reads_overlaps.cbegin(), reads_overlaps.cend(),
                      detail::CmpReadsOvlpsPairsByReadId)) {
    throw std::runtime_error(
        "[camel::CalculateCoverage] input is exptected to be sorted");
  }

  if (std::filesystem::exists(pile_storage_dir)) {
    std::filesystem::remove_all(pile_storage_dir);
  }

  std::filesystem::create_directory(pile_storage_dir);
  auto timer = biosoup::Timer();

  auto const update_coverage =
      [&reads_overlaps](Pile& pile, biosoup::Overlap const& ovlp) -> void {
    auto [query_substr, target_substr] =
        detail::OverlapStrings(reads_overlaps, ovlp);
    auto const edlib_res_align =
        detail::AlignStrings(query_substr, target_substr);

    query_substr.clear();
    query_substr.shrink_to_fit();

    auto query_pos = ovlp.lhs_begin;
    auto target_pos = 0U;

    for (auto i = 0U; i < edlib_res_align.alignmentLength; ++i) {
      switch (edlib_res_align.alignment[i]) {
        case 0:
        case 3: {  // mismatch
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

  auto coverage_futures = std::vector<std::future<Pile>>();

  auto const calc_coverage =
      [&reads_overlaps,
       &update_coverage](ReadOverlapsPair const& read_overlaps) -> Pile {
    auto dst =
        Pile{.id = read_overlaps.read->id,
             .seq_name = read_overlaps.read->name,
             .covgs = std::vector<Coverage>(read_overlaps.read->inflated_len)};

    auto const& local_ovlps = read_overlaps.overlaps;
    for (auto const& target_ovlp : local_ovlps) {
      auto const query_ovlp = detail::ReverseOverlap(target_ovlp);
      update_coverage(dst, query_ovlp);
    }

    for (auto const& [other_read, remote_ovlps] : reads_overlaps) {
      if (read_overlaps.read->id != other_read->id) {
        auto iter = std::lower_bound(
            remote_ovlps.cbegin(), remote_ovlps.cend(), read_overlaps.read->id,
            [](biosoup::Overlap const& ovlp, std::uint32_t const query_id) {
              return ovlp.lhs_id < query_id;
            });

        for (; iter != remote_ovlps.cend() &&
               iter->lhs_id == read_overlaps.read->id;
             ++iter) {
          update_coverage(dst, *iter);
        }
      }
    }

    return dst;
  };

  auto const calc_coverage_async =
      [&thread_pool, &calc_coverage](
          std::reference_wrapper<ReadOverlapsPair const> read_overlaps)
      -> std::future<Pile> {
    return thread_pool->Submit(calc_coverage, read_overlaps);
  };

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

  {
    auto batch_id = 0U;
    auto active_piles = std::vector<Pile>();
    for (auto align_first = reads_overlaps.cbegin();
         align_first != reads_overlaps.cend(); ++batch_id) {
      timer.Start();
      auto const align_last = detail::FindBatchLast(
          align_first, reads_overlaps.cend(), detail::kAlignBatchCap);

      auto const batch_sz = std::distance(align_first, align_last);
      coverage_futures.reserve(batch_sz);
      active_piles.reserve(batch_sz);

      while (!serialize_futures.empty() &&
             is_ser_future_ready(serialize_futures.front())) {
        pop_and_report();
      }

      std::transform(align_first, align_last,
                     std::back_inserter(coverage_futures), calc_coverage_async);

      fmt::print(
          stderr,
          "[camel::CalculateCoverage] calculating coverage for {} reads\n",
          batch_sz);

      std::transform(coverage_futures.begin(), coverage_futures.end(),
                     std::back_inserter(active_piles),
                     std::mem_fn(&std::future<Pile>::get));

      serialize_futures.emplace_back(thread_pool->Submit(
          [&pile_storage_dir](std::vector<Pile> batch_piles,
                              std::size_t batch_id)
              -> std::tuple<std::filesystem::path, std::size_t, std::size_t,
                            double> {
            auto ser_timer = biosoup::Timer();

            ser_timer.Start();
            auto const dst_file = SerializePileBatch(
                batch_piles.cbegin(), batch_piles.cend(), pile_storage_dir,
                fmt::format("pile_batch_{:04d}", batch_id));

            auto const uncompressed_bytes = std::transform_reduce(
                std::make_move_iterator(batch_piles.begin()),
                std::make_move_iterator(batch_piles.end()), 0UL,
                std::plus<std::size_t>(), [](Pile&& pile) -> std::size_t {
                  return sizeof(Pile::id) + pile.seq_name.size() +
                         pile.covgs.size() * sizeof(Coverage);
                });

            auto const compressed_bytes = std::filesystem::file_size(dst_file);

            return std::make_tuple(std::move(dst_file), uncompressed_bytes,
                                   compressed_bytes, ser_timer.Stop());
          },
          std::move(active_piles), batch_id));

      fmt::print(stderr,
                 "[camel::CalculateCoverage]({:12.3f}) finished coverage "
                 "batch; {:3.2f}\% reads covered\n",
                 timer.Stop(),
                 100. * std::distance(reads_overlaps.cbegin(), align_last) /
                     reads_overlaps.size());

      coverage_futures.clear();
      align_first = align_last;
    }
  }

  timer.Start();
  while (!serialize_futures.empty()) {
    pop_and_report();
  }
  timer.Stop();

  fmt::print(stderr, "[camel::CalculateCoverage]({:12.3f})\n",
             timer.elapsed_time());
}

CAMEL_EXPORT auto CallBases(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::vector<ReadOverlapsPair> const& reads_overlaps)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
  if (!std::is_sorted(reads_overlaps.cbegin(), reads_overlaps.cend(),
                      detail::CmpReadsOvlpsPairsByReadId)) {
    throw std::runtime_error("[camel::CallBases] int is expected to be sorted");
  }

  auto dst = std::vector<std::unique_ptr<biosoup::NucleicAcid>>();
  dst.reserve(reads_overlaps.size());

  auto const update_coverage =
      [&reads_overlaps](
          std::vector<Coverage>& base_covgs,
          tsl::robin_map<std::uint32_t, std::vector<detail::InsSignal>>&
              insertions,
          biosoup::Overlap const& ovlp, std::uint16_t ins_freq) -> void {
    auto [query_substr, target_substr] =
        detail::OverlapStrings(reads_overlaps, ovlp);

    auto const edlib_res_align =
        detail::AlignStrings(query_substr, target_substr);

    auto query_pos = ovlp.lhs_begin;
    auto target_pos = 0U;

    for (auto i = 0U; i < edlib_res_align.alignmentLength; ++i) {
      switch (edlib_res_align.alignment[i]) {
        case 0:
        case 3: {
          /* clang-format off */
          switch (target_substr[target_pos]) {
            case 'A': ++base_covgs[query_pos].a;
            case 'C': ++base_covgs[query_pos].c;
            case 'G': ++base_covgs[query_pos].g;
            case 'T': ++base_covgs[query_pos].t;
            default: break;
          }
          /* clang-format on */

          ++query_pos;
          ++target_pos;

          break;
        }
        case 1: {  // insertion on target
          ++base_covgs[query_pos].del;
          ++query_pos;

          break;
        }
        case 2: {  // insertion on query
          if (base_covgs[query_pos].ins < 10) {
            ++base_covgs[query_pos].ins;
            ++target_pos;
          } else {
            auto& ins_vec = insertions[query_pos];
            if (ins_vec.empty()) {
              ins_vec.reserve(4U);
            }

            for (auto j = 0U; i < edlib_res_align.alignmentLength &&
                              edlib_res_align.alignment[i] == 2;
                 ++j, ++i) {
              if (j >= ins_vec.size()) {
                ins_vec.emplace_back();
              }

              /* clang-format off */
              switch (target_substr[target_pos]) {
                case 'A': ++ins_vec[j].a; break;
                case 'C': ++ins_vec[j].c; break;
                case 'G': ++ins_vec[j].g; break;
                case 'T': ++ins_vec[j].t; break;
                default: break;
              }
              /* clang-format on */

              ++base_covgs[query_pos].ins;
              ++target_pos;
            }

            --i;
          }

          break;
        }

        default: {
          break;
        }
      }
    }
  };

  auto call_bases = [&reads_overlaps,
                     &update_coverage](ReadOverlapsPair const& read_overlaps)
      -> std::unique_ptr<biosoup::NucleicAcid> {
    auto dst = std::unique_ptr<biosoup::NucleicAcid>();

    auto base_covgs = std::vector<Coverage>(read_overlaps.read->inflated_len);
    auto insertions =
        tsl::robin_map<std::uint32_t, std::vector<detail::InsSignal>>();

    auto const& local_ovlsp = read_overlaps.overlaps;
    for (auto const& target_ovlp : local_ovlsp) {
      auto const query_ovlp = detail::ReverseOverlap(target_ovlp);
      // TODO: dynamically determine ins freq threshold
      update_coverage(base_covgs, insertions, query_ovlp, 10U);
    }

    for (auto const& [other_read, remote_ovlps] : reads_overlaps) {
      if (read_overlaps.read->id != other_read->id) {
        auto iter = std::lower_bound(
            remote_ovlps.cbegin(), remote_ovlps.cend(), read_overlaps.read->id,
            [](biosoup::Overlap const& ovlp, std::uint32_t const query_id) {
              return ovlp.lhs_id < query_id;
            });

        for (; iter != remote_ovlps.cend() &&
               iter->lhs_id == read_overlaps.read->id;
             ++iter) {
          update_coverage(base_covgs, insertions, *iter, 10U);
        }
      }
    }

    // call bases
    {
      auto build_buff = std::string();
      build_buff.reserve(read_overlaps.read->inflated_len * 1.1);

      for (auto pos = 0U; pos < base_covgs.size(); ++pos) {
        auto base_val = char(0);
        auto const kBaseCovg = base_covgs[pos].a + base_covgs[pos].c +
                               base_covgs[pos].g + base_covgs[pos].t;

        if (base_covgs[pos].del >= kBaseCovg) {
          base_val = 0;
        } else if (2 * base_covgs[pos].a > kBaseCovg) {
          base_val = 'A';
        } else if (2 * base_covgs[pos].c > kBaseCovg) {
          base_val = 'C';
        } else if (2 * base_covgs[pos].g > kBaseCovg) {
          base_val = 'G';
        } else if (2 * base_covgs[pos].t > kBaseCovg) {
          base_val = 'T';
        } else {
          base_val = biosoup::kNucleotideDecoder[read_overlaps.read->Code(pos)];
        }

        if (base_val > 0) {
          build_buff.push_back(base_val);
        }

        if (1.2 * base_covgs[pos].ins >= kBaseCovg) {
          for (auto const it : insertions[pos]) {
            auto max_val = std::max({it.a, it.c, it.g, it.t});
            if (max_val == it.a) {
              build_buff.push_back('A');
              continue;
            }
            if (max_val == it.c) {
              build_buff.push_back('C');
              continue;
            }
            if (max_val == it.g) {
              build_buff.push_back('G');
              continue;
            }
            if (max_val == it.t) {
              build_buff.push_back('T');
              continue;
            }
          }
        }
      }

      dst = std::make_unique<biosoup::NucleicAcid>(read_overlaps.read->name,
                                                   build_buff);
    }

    return dst;
  };

  auto base_call_futures =
      std::vector<std::future<std::unique_ptr<biosoup::NucleicAcid>>>();

  auto const call_bases_async =
      [&thread_pool, &call_bases](
          std::reference_wrapper<ReadOverlapsPair const> read_overlaps)
      -> std::future<std::unique_ptr<biosoup::NucleicAcid>> {
    return thread_pool->Submit(call_bases, read_overlaps);
  };

  auto timer = biosoup::Timer();
  for (auto batch_first = reads_overlaps.cbegin();
       batch_first != reads_overlaps.cend();) {
    timer.Start();
    auto const batch_last = detail::FindBatchLast(
        batch_first, reads_overlaps.cend(), detail::kAlignBatchCap);

    auto const batch_sz = std::distance(batch_first, batch_last);
    base_call_futures.reserve(batch_sz);

    std::transform(batch_first, batch_last,
                   std::back_inserter(base_call_futures), call_bases_async);

    fmt::print(stderr,
               "[camel::CallBases] estimateing correct bases for {} reads\n",
               batch_sz);

    std::transform(
        base_call_futures.begin(), base_call_futures.end(),
        std::back_inserter(dst),
        std::mem_fn(&std::future<std::unique_ptr<biosoup::NucleicAcid>>::get));

    fmt::print(stderr,
               "[camel::CallBases]({:12.3f}) finisehd base call batch; "
               "{:3.2f}\% reads covered\n",
               timer.Stop(),
               100. * std::distance(reads_overlaps.cbegin(), batch_last) /
                   reads_overlaps.size());

    base_call_futures.clear();
    batch_first = batch_last;
  }

  fmt::print(stderr, "[camel::CallBases]({:12.3f})\n", timer.elapsed_time());

  return dst;
}

}  // namespace camel
