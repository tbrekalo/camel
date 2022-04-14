#include "camel/io.h"

#include <array>
#include <fstream>
#include <functional>
#include <iterator>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string_view>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "detail/serialization.h"
#include "fmt/core.h"

namespace camel {

namespace detail {

static constexpr std::size_t kDefaultPileStorageFileSz = 1UL << 32UL;  // ~4gb

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

}  // namespace detail

auto LoadSequences(std::filesystem::path const& path)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
  auto dst = std::vector<std::unique_ptr<biosoup::NucleicAcid>>();

  auto parser = detail::CreateParser(path);
  auto local_reads = parser->Parse(std::numeric_limits<std::uint64_t>::max());

  std::move(local_reads.begin(), local_reads.end(), std::back_inserter(dst));
  dst.shrink_to_fit();

  return dst;
}

auto LoadSequences(std::shared_ptr<thread_pool::ThreadPool> thread_pool,
                   std::vector<std::filesystem::path> const& paths)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
  auto parse_futures = std::vector<
      std::future<std::vector<std::unique_ptr<biosoup::NucleicAcid>>>>();
  parse_futures.reserve(paths.size());

  std::transform(
      paths.cbegin(), paths.cend(), std::back_inserter(parse_futures),
      [&thread_pool](std::filesystem::path const& path) {
        return thread_pool->Submit(
            [](std::filesystem::path const& p) { return LoadSequences(p); },
            path);
      });

  auto dst = std::vector<std::unique_ptr<biosoup::NucleicAcid>>();
  for (auto& it : parse_futures) {
    auto&& vec = it.get();
    std::move(vec.begin(), vec.end(), std::back_inserter(dst));
  }

  dst.shrink_to_fit();
  return dst;
}

template <class Buff>
static auto Store(Buff& buff, Pile const& pile) -> void {
  buff(pile.id, pile.seq_name, pile.covgs);
}

template <class Buff>
static auto Load(Buff& buff, Pile& pile) -> void {
  buff(pile.id, pile.seq_name, pile.covgs);
}

CAMEL_EXPORT auto SerializePileBatch(std::vector<Pile>::const_iterator first,
                                     std::vector<Pile>::const_iterator last,
                                     std::filesystem::path const& dst_dir,
                                     std::string const& batch_name)
    -> std::filesystem::path {
  auto de_sz = std::size_t(0);

  auto binary_out = detail::BinaryOutBuffer();
  binary_out(static_cast<std::size_t>(std::distance(first, last)));
  for (auto it = first; it != last; ++it) {
    de_sz += sizeof(it->id) + it->seq_name.size() +
             it->covgs.size() * sizeof(Coverage);
    binary_out(*it);
  }

  auto const dst_path = dst_dir / (batch_name + ".camel");

  detail::CompressToFile(binary_out.Bytes(), dst_path);
  return dst_path;
}

auto SerializePiles(std::shared_ptr<thread_pool::ThreadPool> thread_pool,
                    std::vector<Pile> const& piles,
                    std::filesystem::path const& dst_dir) -> void {
  SerializePiles(std::move(thread_pool), piles, dst_dir,
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
    auto const pile_file_name = fmt::format("pile_dump_{:04d}", file_id);
    SerializePileBatch(first, last, dst_dir, pile_file_name);
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

          auto binary_in =
              detail::BinaryInBuffer(detail::DecompressFromFile(dir_entry));

          auto sz = std::size_t();
          binary_in(sz);

          dst.resize(sz);
          for (auto i = 0UL; i < sz; ++i) {
            binary_in(dst[i]);
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

  std::sort(
      dst.begin(), dst.end(),
      [](Pile const& lhs, Pile const& rhs) -> bool { return lhs.id < rhs.id; });

  dst.shrink_to_fit();
  return dst;
}

}  // namespace camel
