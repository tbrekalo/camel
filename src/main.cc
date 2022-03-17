#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

#include "bioparser/fastq_parser.hpp"
#include "biosoup/nucleic_acid.hpp"
#include "cxxopts.hpp"
#include "fmt/core.h"
#include "thread_pool/thread_pool.hpp"

int main(int argc, char** argv) {
  auto options =
      cxxopts::Options("camel", "Camel is haplotype aware detection tool");
  options.add_options()("t,threads", "number of threads avalable for execution",
                        cxxopts::value<std::uint32_t>())(
      "reads", "input fastq reads", cxxopts::value<std::vector<std::string>>());

  options.parse_positional("reads");
  auto const result = options.parse(argc, argv);

  auto thread_pool = std::make_shared<thread_pool::ThreadPool>(
      result["threads"].as<std::uint32_t>());

  auto const fastq_reads = result["reads"].as<std::vector<std::string>>();

  return EXIT_SUCCESS;
}
