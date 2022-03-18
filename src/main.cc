#include <algorithm>
#include <array>
#include <cstdlib>
#include <functional>
#include <iterator>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "bioparser/fastq_parser.hpp"
#include "biosoup/nucleic_acid.hpp"
#include "cxxopts.hpp"
#include "fmt/core.h"
#include "thread_pool/thread_pool.hpp"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects(0);

static auto IsFastqFilename(std::string_view const filename) -> bool {
  auto const is_suffix = [](std::string_view const suffix,
                            std::string_view const target) -> bool {
    return suffix.length() <= target.length()
               ? suffix == target.substr(target.length() - suffix.length())
               : false;
  };

  static constexpr auto kFastqSuffixes =
      std::array<char const*, 4>{".fastq", ".fq", "fastq.gz", "fq.gz"};

  using namespace std::placeholders;

  return std::any_of(kFastqSuffixes.cbegin(), kFastqSuffixes.end(),
                     std::bind(is_suffix, _1, filename));
}

static auto CreateFastqParser(std::string const& path)
    -> std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>> {
  return bioparser::Parser<biosoup::NucleicAcid>::Create<
      bioparser::FastqParser>(path);
}

auto main(int argc, char** argv) -> int {
  auto options =
      cxxopts::Options("camel", "Camel is haplotype aware detection tool");
  options.add_options()("t,threads", "number of threads avalable for execution",
                        cxxopts::value<std::uint32_t>())(
      "reads", "input fastq reads", cxxopts::value<std::vector<std::string>>());

  options.parse_positional("reads");
  auto const result = options.parse(argc, argv);

  auto thread_pool = std::make_shared<thread_pool::ThreadPool>(
      result["threads"].as<std::uint32_t>());

  auto reads = std::vector<std::unique_ptr<biosoup::NucleicAcid>>();

  {
    auto fastq_read_path_str = result["reads"].as<std::vector<std::string>>();
    if (!std::all_of(fastq_read_path_str.cbegin(), fastq_read_path_str.cend(),
                     IsFastqFilename)) {
      auto const invalid_begin =
          std::partition(fastq_read_path_str.begin(), fastq_read_path_str.end(),
                         IsFastqFilename);

      std::for_each(invalid_begin, fastq_read_path_str.end(),
                    [](std::string_view filename) -> void {
                      fmt::print(stderr, "[camel] invalid file name {}\n",
                                 filename);
                    });

      return EXIT_FAILURE;
    }
    auto parsers =
        std::vector<std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>>>();
    parsers.reserve(fastq_read_path_str.size());

    std::transform(fastq_read_path_str.cbegin(), fastq_read_path_str.cend(),
                   std::back_inserter(parsers), CreateFastqParser);

    for (auto& parser : parsers) {
      auto local_reads =
          parser->Parse(std::numeric_limits<std::uint64_t>::max());

      std::move(local_reads.begin(), local_reads.end(),
                std::back_inserter(reads));
    }
  }

  fmt::print(stderr, "[camel] loaded {} reads\n", reads.size());

  return EXIT_SUCCESS;
}
