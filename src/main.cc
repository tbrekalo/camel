#include <algorithm>
#include <array>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

#include "camel/io.h"
#include "cxxopts.hpp"
#include "fmt/core.h"
#include "thread_pool/thread_pool.hpp"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects(0);

auto main(int argc, char** argv) -> int {
  auto options =
      cxxopts::Options("camel", "Camel is haplotype aware detection tool");
  options.add_options()("t,threads", "number of threads avalable for execution",
                        cxxopts::value<std::uint32_t>())(
      "paths", "input fastq reads", cxxopts::value<std::vector<std::string>>());

  options.parse_positional("paths");
  auto const result = options.parse(argc, argv);

  auto thread_pool = std::make_shared<thread_pool::ThreadPool>(
      result["threads"].as<std::uint32_t>());

  auto paths = std::vector<std::filesystem::path>();
  auto paths_strs = result["paths"].as<std::vector<std::string>>();
  std::transform(std::make_move_iterator(paths_strs.begin()),
                 std::make_move_iterator(paths_strs.end()),
                 std::back_inserter(paths),
                 [](std::string&& path_str) -> std::filesystem::path {
                   return std::filesystem::path(std::move(path_str));
                 });

  auto reads = camel::LoadSequences(thread_pool, paths);
  fmt::print(stderr, "[camel] loaded {} reads\n", reads.size());

  return EXIT_SUCCESS;
}
