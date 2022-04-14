#include <algorithm>
#include <array>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

#include "biosoup/timer.hpp"
#include "camel/coverage.h"
#include "camel/io.h"
#include "camel/mapping.h"
#include "cxxopts.hpp"
#include "fmt/core.h"
#include "thread_pool/thread_pool.hpp"

#include "mimalloc-new-delete.h"
#include "mimalloc-override.h"
#include "mimalloc.h"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects(0);

auto main(int argc, char** argv) -> int {
  auto options =
      cxxopts::Options("camel", "Camel is haplotype aware detection tool");
  options.add_options()("t,threads", "number of threads avalable for execution",
                        cxxopts::value<std::uint32_t>())(
      "s,serialization_dst", "destination folder for pile serialization",
      cxxopts::value<std::string>())(
      "paths", "input fastq reads", cxxopts::value<std::vector<std::string>>());
  options.parse_positional("paths");

  auto const result = options.parse(argc, argv);

  auto thread_pool = std::make_shared<thread_pool::ThreadPool>(
      result["threads"].as<std::uint32_t>());

  auto ser_dst_path =
      std::filesystem::path(result["serialization_dst"].as<std::string>());

  auto paths = std::vector<std::filesystem::path>();
  auto paths_strs = result["paths"].as<std::vector<std::string>>();
  std::transform(std::make_move_iterator(paths_strs.begin()),
                 std::make_move_iterator(paths_strs.end()),
                 std::back_inserter(paths),
                 [](std::string&& path_str) -> std::filesystem::path {
                   return std::filesystem::path(std::move(path_str));
                 });

  auto timer = biosoup::Timer();

  timer.Start();
  auto const reads = camel::LoadSequences(thread_pool, paths);
  fmt::print(stderr, "[camel]({:12.3f}) loaded {} reads\n", timer.Stop(),
             reads.size());
  timer.Start();

  camel::CalculateCoverage(thread_pool, camel::MapCfg{}, reads, ser_dst_path);
  timer.Stop();

  fmt::print(stderr, "[camel]({:12.3f}) done\n", timer.elapsed_time());

  return EXIT_SUCCESS;
}
