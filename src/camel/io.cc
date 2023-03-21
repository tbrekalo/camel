#include "camel/io.h"

#include <fstream>
#include <unordered_map>

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

namespace camel {

namespace detail {

static constexpr std::size_t kDefaultPileStorageFileSz = 1UL << 32UL;  // 4GiB
static constexpr std::size_t kDefaultSeqStorageFileSz = 1UL << 22UL;   // 4MiB

static constexpr auto kFastaSuffxies =
    std::array<char const*, 4>{".fasta", "fasta.gz", ".fa", ".fa.gz"};

static constexpr auto kFastqSuffixes =
    std::array<char const*, 4>{".fastq", ".fastq.gz", ".fq", ".fq.gz"};

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
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& seqs,
    std::filesystem::path const& dst_file) -> void {
  if (seqs.empty()) {
    return;
  }

  auto ofstrm =
      std::fstream(dst_file, std::ios_base::out | std::ios_base::trunc);
  SerializeSequences(seqs, ofstrm);
}

CAMEL_EXPORT auto SerializeSequences(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& seqs,
    std::ostream& ofstrm) -> void {
  for (auto const& seq : seqs) {
    ofstrm << '>' << seq->name << '\n' << seq->InflateData() << '\n';
  }
}

auto LoadSequences(std::vector<std::filesystem::path> const& paths)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
  auto buff_vec =
      tbb::concurrent_vector<std::unique_ptr<biosoup::NucleicAcid>>();
  auto dst = std::vector<std::unique_ptr<biosoup::NucleicAcid>>();

  tbb::parallel_for_each(
      paths, [&buff_vec](std::filesystem::path const& path) -> void {
        for (auto& seq : LoadSequences(path)) {
          buff_vec.push_back(std::move(seq));
        }
      });

  // reindexing
  tbb::parallel_sort(buff_vec, detail::CmpNucleicAcidByName);
  tbb::parallel_for(
      0U, static_cast<std::uint32_t>(buff_vec.size()),
      [&](std::uint32_t idx) -> void { buff_vec[idx]->id = idx; });

  dst.reserve(buff_vec.size());
  std::move(buff_vec.begin(), buff_vec.end(), std::back_inserter(dst));

  return dst;
}

auto LoadOverlaps(
    std::filesystem::path const& paf_path,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads,
    double const error_threshold, std::size_t const n_overlaps)
    -> std::vector<std::vector<biosoup::Overlap>> {
  auto dst = std::vector<std::vector<biosoup::Overlap>>(reads.size());
  auto name_to_id =
      std::unordered_map<std::string_view, std::uint32_t>(reads.size());

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
    auto parser =
        bioparser::Parser<detail::PafOverlap>::Create<bioparser::PafParser>(
            paf_path.string());

    auto ovlps = parser->Parse(std::numeric_limits<std::uint64_t>::max(), true);

    for (auto& ovlp_ptr : ovlps) {
      auto ovlp = transform_overlap(std::move(ovlp_ptr));
      if (ovlp.lhs_id != ovlp.rhs_id &&
          detail::OverlapError(ovlp) <= error_threshold) {
        dst[ovlp.rhs_id].push_back(ovlp);
      }
    }
  }

  return dst;
}

}  // namespace camel
