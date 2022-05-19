#include <algorithm>
#include <array>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

#include "biosoup/timer.hpp"
#include "camel/correct.h"
#include "camel/coverage.h"
#include "camel/detail/overlap.h"
#include "camel/io.h"
#include "camel/mapping.h"
#include "camel/state.h"
#include "cxxopts.hpp"
#include "fmt/core.h"
#include "fmt/ostream.h"
// #include "jemalloc/jemalloc.h"
#include "thread_pool/thread_pool.hpp"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects(0);

auto main(int argc, char** argv) -> int {
  auto options =
      cxxopts::Options("camel", "Camel is haplotype aware detection tool");

  /* clang-format off */
  options.add_options()
    ("t,threads", "number of threads avalable for execution",
            cxxopts::value<std::uint32_t>())
    ("s,serialization_dst", "destination folder for pile serialization",
            cxxopts::value<std::string>()->default_value("./camel_piles"))
    ("l,log_dst", "destination for for logging information",
            cxxopts::value<std::string>()->default_value("./camel_log"))
    ("o,out", "output destination folder",
      cxxopts::value<std::string>()->default_value("./camel_out"))
    ("paths", "input fastq reads", 
            cxxopts::value<std::vector<std::string>>());
  options.parse_positional("paths");
  /* clang-format on */

  auto const result = options.parse(argc, argv);

  auto state =
      camel::State{.thread_pool = std::make_shared<thread_pool::ThreadPool>(
                       result["threads"].as<std::uint32_t>()),
                   .log_path = result["log_dst"].as<std::string>()};

  auto const ser_dst_path =
      std::filesystem::path(result["serialization_dst"].as<std::string>());
  auto const out_path = std::filesystem::path(result["out"].as<std::string>());

  if (std::filesystem::exists(state.log_path)) {
    std::filesystem::remove_all(state.log_path);
  }

  std::filesystem::create_directory(state.log_path);

  if (std::filesystem::exists(out_path)) {
    std::filesystem::remove_all(out_path);
  }

  std::filesystem::create_directory(out_path);

  auto paths = std::vector<std::filesystem::path>();

  {
    auto paths_strs = result["paths"].as<std::vector<std::string>>();
    for (auto const& it : paths_strs) {
      auto it_path = std::filesystem::path(std::move(it));
      if (std::filesystem::is_regular_file(it_path)) {
        paths.emplace_back(std::move(it_path));
      } else if (std::filesystem::is_directory(it_path)) {
        for (auto dir_entry : std::filesystem::recursive_directory_iterator(
                 std::move(it_path))) {
          if (std::filesystem::is_regular_file(dir_entry)) {
            paths.emplace_back(std::move(dir_entry));
          }
        }
      }
    }
  }

  decltype(paths)(std::make_move_iterator(paths.begin()),
                  std::make_move_iterator(paths.end()))
      .swap(paths);

  auto timer = biosoup::Timer();

  {
    timer.Start();
    auto reads = camel::LoadSequences(state, paths);
    fmt::print(stderr, "[camel]({:12.3f}) loaded {} reads\n", timer.Stop(),
               reads.size());

    timer.Start();
    auto corrected_and_annoted =
        camel::SnpErrorCorrect(state, camel::PolishConfig{}, std::move(reads));
    timer.Stop();

    auto dump = std::vector<std::unique_ptr<biosoup::NucleicAcid>>();

    dump.reserve(corrected_and_annoted.size());
    std::transform(std::make_move_iterator(corrected_and_annoted.begin()),
                   std::make_move_iterator(corrected_and_annoted.end()),
                   std::back_inserter(dump),
                   [](camel::AnnotatedRead ar)
                       -> std::unique_ptr<biosoup::NucleicAcid> {
                     return std::move(ar.read);
                   });

    camel::StoreSequences(state, dump, out_path, 1U << 28U);
  }

  fmt::print(stderr, "[camel]({:12.3f}) done\n", timer.elapsed_time());

  return EXIT_SUCCESS;
}
