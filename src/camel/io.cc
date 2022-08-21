#include "camel/io.h"

#include <array>
#include <fstream>
#include <functional>
#include <future>
#include <iterator>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string_view>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "bioparser/paf_parser.hpp"
#include "detail/overlap.h"
#include "fmt/core.h"
#include "tbb/concurrent_vector.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_for_each.h"
#include "tbb/parallel_sort.h"
#include "tbb/task_group.h"
#include "tsl/robin_map.h"

namespace camel {

namespace detail {

static constexpr std::size_t kDefaultPileStorageFileSz = 1UL << 32UL;  // 4GiB
static constexpr std::size_t kDefaultSeqStorageFileSz = 1UL << 22UL;   // 4MiB

static constexpr auto kFastaSuffxies =
    std::array<char const*, 4>{".fasta", "fasta.gz", ".fa", ".fa.gz"};

static constexpr auto kFastqSuffixes =
    std::array<char const*, 4>{".fastq", ".fastq.gz", ".fa", ".fa.gz"};

static auto IsSuffixFor(std::string_view const suffix,
                        std::string_view const query) -> bool {
  return suffix.length() <= query.length()
             ? suffix == query.substr(query.length() - suffix.length())
             : false;
}

static auto CreateParser(std::filesystem::path const& path)
    -> std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>> {
  using namespace std::placeholders;
  if (std::filesystem::exists(path)) {
    if (std::any_of(kFastaSuffxies.cbegin(), kFastaSuffxies.cend(),
                    std::bind(IsSuffixFor, _1, path.c_str()))) {
      return bioparser::Parser<biosoup::NucleicAcid>::Create<
          bioparser::FastaParser>(path.c_str());
    } else if (std::any_of(kFastqSuffixes.cbegin(), kFastqSuffixes.cend(),
                           std::bind(IsSuffixFor, _1, path.c_str()))) {
      return bioparser::Parser<biosoup::NucleicAcid>::Create<
          bioparser::FastqParser>(path.c_str());
    }
  }

  throw std::invalid_argument(
      "[camel::detail::CreateParser] invalid file path: " + path.string());
}

static auto CmpNucleicAcidByName(
    std::unique_ptr<biosoup::NucleicAcid> const& lhs,
    std::unique_ptr<biosoup::NucleicAcid> const& rhs) -> bool {
  return lhs->name < rhs->name;
}

/* clang-format off */
struct PafOverlap {

  PafOverlap() = default;

  PafOverlap(PafOverlap const&) = default;
  PafOverlap& operator=(PafOverlap const&) = default;

  PafOverlap(PafOverlap&&) = default;
  PafOverlap& operator=(PafOverlap&&) = default;

  PafOverlap(
    char const* query_name,
    std::uint32_t query_name_len,

    std::uint32_t query_len,
    std::uint32_t query_start,
    std::uint32_t query_end,

    char relative_strand,

    char const* target_name,
    std::uint32_t target_name_len,

    std::uint32_t target_len,
    std::uint32_t target_start,
    std::uint32_t target_end,
    
    std::uint32_t num_matches,
    std::uint32_t alignment_len,
    std::uint32_t mapping_quality) :
      query_name(query_name, query_name_len),
      query_len(query_len),
      query_start(query_start),
      query_end(query_end),

      relative_strand(relative_strand),

      target_name(target_name, target_name_len),
      target_len(target_len),
      target_start(target_start),
      target_end(target_end),

      num_matches(num_matches),
      alignment_len(alignment_len),
      mapping_quality(mapping_quality) {}
  
  std::string query_name;
  std::uint32_t query_len;
  std::uint32_t query_start;
  std::uint32_t query_end;

  char relative_strand;

  std::string target_name;
  std::uint32_t target_len;
  std::uint32_t target_start;
  std::uint32_t target_end;

  std::uint32_t num_matches;
  std::uint32_t alignment_len;
  std::uint32_t mapping_quality;
};
/* clang-format on */

}  // namespace detail

auto LoadSequences(std::filesystem::path const& path)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
  auto parser = detail::CreateParser(path);
  auto dst = parser->Parse(std::numeric_limits<std::uint64_t>::max());

  std::sort(dst.begin(), dst.end(), detail::CmpNucleicAcidByName);

  // reindexing
  for (auto idx = 0U; idx < dst.size(); ++idx) {
    dst[idx]->id = idx;
  }

  decltype(dst)(std::make_move_iterator(dst.begin()),
                std::make_move_iterator(dst.end()))
      .swap(dst);

  return dst;
}

auto StoreSequences(
    tbb::task_arena& task_arena,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& seqs,
    std::filesystem::path const& dst_folder) -> void {
  StoreSequences(task_arena, seqs, dst_folder,
                 detail::kDefaultSeqStorageFileSz);
}

auto StoreSequences(
    tbb::task_arena& task_arena,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& seqs,
    std::filesystem::path const& dst_folder, std::uint64_t dst_file_cap)
    -> void {
  if (seqs.empty()) {
    return;
  }

  auto const find_batch_last =
      [](std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator
             first,
         std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator
             last,
         std::uint64_t const batch_cap)
      -> std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator {
    for (auto batch_sz = 0UL; batch_sz < batch_cap && first != last; ++first) {
      batch_sz += 2UL * ((*first)->inflated_len);
    }

    return first;
  };

  auto const is_fasta = seqs.front()->block_quality.empty();

  auto const store_fasta_impl =
      +[](std::ostream& ostrm,
          std::unique_ptr<biosoup::NucleicAcid> const& seq) -> void {
    ostrm << '>' << seq->name << '\n' << seq->InflateData() << '\n';
  };

  auto const store_fastq_impl =
      +[](std::ostream& ostrm,
          std::unique_ptr<biosoup::NucleicAcid> const& seq) -> void {
    ostrm << '@' << seq->name << '\n'
          << seq->InflateData() << '\n'
          << "+\n"
          << seq->InflateQuality() << '\n';
  };

  using StoreFnSig =
      void(std::ostream&, std::unique_ptr<biosoup::NucleicAcid> const&);
  using StoreFnPtr = std::add_const_t<std::add_pointer_t<StoreFnSig>>;

  static_assert(std::is_same_v<decltype(store_fasta_impl), StoreFnPtr>);
  static_assert(std::is_same_v<decltype(store_fastq_impl), StoreFnPtr>);

  auto const store_fn_impl = is_fasta ? store_fasta_impl : store_fastq_impl;

  auto const store_fn =
      [impl = store_fn_impl](
          std::fstream& ofstrm,
          std::unique_ptr<biosoup::NucleicAcid> const& seq) -> void {
    impl(ofstrm, seq);
  };

  {
    task_arena.execute([&]() -> void {
      auto batch_id = 0U;
      auto task_group = tbb::task_group();
      for (auto first = seqs.cbegin(); first != seqs.cend(); ++batch_id) {
        auto const last = find_batch_last(first, seqs.cend(), dst_file_cap);
        task_group.run([&dst_folder, is_fasta, store_fn, first, last,
                        batch_id]() -> void {
          auto local_first = first;
          auto local_last = last;
          auto const dst_file_path =
              dst_folder / (fmt::format("corrected_batch_{:04d}", batch_id) +
                            (is_fasta ? ".fa" : ".fq"));

          auto ofstrm =
              std::fstream(dst_file_path, std::ios::out | std::ios::trunc);

          for (; local_first != local_last; ++local_first) {
            store_fn(ofstrm, *first);
          }
        });

        first = last;
      }

      task_group.wait();
    });
  }
}

auto LoadSequences(tbb::task_arena& task_arena,
                   std::vector<std::filesystem::path> const& paths)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
  auto buff_vec =
      tbb::concurrent_vector<std::unique_ptr<biosoup::NucleicAcid>>();
  auto dst = std::vector<std::unique_ptr<biosoup::NucleicAcid>>();

  task_arena.execute([&buff_vec, &dst, &paths]() -> void {
    tbb::parallel_for_each(
        paths, [&buff_vec](std::filesystem::path const& path) -> void {
          for (auto& seq : LoadSequences(path)) {
            buff_vec.push_back(std::move(seq));
          }
        });

    // reindexing
    tbb::parallel_sort(buff_vec);
    tbb::parallel_for(
        0U, static_cast<std::uint32_t>(buff_vec.size()),
        [&](std::uint32_t idx) -> void { buff_vec[idx]->id = idx; });
  });

  dst.reserve(buff_vec.size());
  std::move(buff_vec.begin(), buff_vec.end(), std::back_inserter(dst));

  return dst;
}

auto LoadOverlaps(
    std::filesystem::path const& paf_path,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads)
    -> std::vector<std::vector<biosoup::Overlap>> {
  auto dst = std::vector<std::vector<biosoup::Overlap>>(reads.size());
  auto name_to_id =
      tsl::robin_map<std::string_view, std::uint32_t>(reads.size());

  for (auto const& read : reads) {
    name_to_id[read->name] = read->id;
  }

  auto const transform_overlap =
      [&name_to_id](
          std::unique_ptr<detail::PafOverlap> paf_ovlp) -> biosoup::Overlap {
    return biosoup::Overlap(name_to_id.at(paf_ovlp->query_name),
                            paf_ovlp->query_start, paf_ovlp->query_end,

                            name_to_id.at(paf_ovlp->target_name),
                            paf_ovlp->target_start, paf_ovlp->target_end,
                            paf_ovlp->num_matches,
                            paf_ovlp->relative_strand == '+');
  };

  {
    auto best_ovlps = std::vector<biosoup::Overlap>(reads.size());
    {
      auto parser =
          bioparser::Parser<detail::PafOverlap>::Create<bioparser::PafParser>(
              paf_path.string());

      auto ovlps =
          parser->Parse(std::numeric_limits<std::uint64_t>::max(), true);

      for (auto& ovlp_ptr : ovlps) {
        auto ovlp = transform_overlap(std::move(ovlp_ptr));
        if (ovlp.lhs_id != ovlp.rhs_id && detail::OverlapLength(ovlp) > 1280U &&
            detail::OverlapError(ovlp) < 0.2 &&
            detail::OverlapLength(best_ovlps[ovlp.lhs_id]) <
                detail::OverlapLength(ovlp)) {
          best_ovlps[ovlp.lhs_id] = ovlp;
        }
      }
    }

    for (auto const& ovlp : best_ovlps) {
      if (detail::OverlapLength(ovlp) > 0U) {
        auto const rev_ovlp = detail::ReverseOverlap(ovlp);
        dst[ovlp.rhs_id].push_back(detail::ReverseOverlap(ovlp));
      }
    }
  }

  return dst;
}

}  // namespace camel
