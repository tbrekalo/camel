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
    State& state,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& seqs,
    std::filesystem::path const& dst_folder) -> void {
  StoreSequences(state, seqs, dst_folder, detail::kDefaultSeqStorageFileSz);
}

auto StoreSequences(
    State& state,
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

  auto ser_futures = std::vector<std::future<void>>();

  {
    auto batch_id = 0U;
    for (auto first = seqs.cbegin(); first != seqs.cend(); ++batch_id) {
      auto const last = find_batch_last(first, seqs.cend(), dst_file_cap);
      ser_futures.emplace_back(state.thread_pool->Submit(
          [&dst_folder, is_fasta, store_fn](
              std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator
                  first,
              std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator
                  last,
              std::uint32_t batch_id) -> void {
            auto const dst_file_path =
                dst_folder / (fmt::format("corrected_batch_{:04d}", batch_id) +
                              (is_fasta ? ".fa" : ".fq"));

            auto ofstrm =
                std::fstream(dst_file_path, std::ios::out | std::ios::trunc);

            for (; first != last; ++first) {
              store_fn(ofstrm, *first);
            }
          },
          first, last, batch_id));

      first = last;
    }
  }

  std::for_each(ser_futures.begin(), ser_futures.end(),
                std::mem_fn(&std::future<void>::get));
}

auto LoadSequences(State& state,
                   std::vector<std::filesystem::path> const& paths)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
  auto parse_futures = std::vector<
      std::future<std::vector<std::unique_ptr<biosoup::NucleicAcid>>>>();
  parse_futures.reserve(paths.size());

  std::transform(
      paths.cbegin(), paths.cend(), std::back_inserter(parse_futures),
      [&thread_pool = state.thread_pool](std::filesystem::path const& path) {
        return thread_pool->Submit(
            [](std::filesystem::path const& p) { return LoadSequences(p); },
            path);
      });

  auto dst = std::vector<std::unique_ptr<biosoup::NucleicAcid>>();
  for (auto& it : parse_futures) {
    auto&& vec = it.get();
    std::move(vec.begin(), vec.end(), std::back_inserter(dst));
  }

  std::sort(dst.begin(), dst.end(),
            [](auto const& lhs, auto const& rhs) -> bool {
              return lhs->id < rhs->id;
            });

  // reindexing
  for (auto idx = 0U; idx < dst.size(); ++idx) {
    dst[idx]->id = idx;
  }

  decltype(dst)(std::make_move_iterator(dst.begin()),
                std::make_move_iterator(dst.end()))
      .swap(dst);

  return dst;
}

}  // namespace camel
