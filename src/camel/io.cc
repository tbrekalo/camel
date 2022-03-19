#include "camel/io.h"

#include <array>
#include <functional>
#include <iterator>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string_view>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"

namespace camel {

namespace detail {

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
      [tp = thread_pool](std::filesystem::path const& path) {
        return tp->Submit(
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

}  // namespace camel
