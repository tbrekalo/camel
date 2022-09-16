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
#include "camel/io.h"
#include "cxxopts.hpp"
#include "fmt/core.h"
#include "tbb/task_arena.h"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects(0);

auto main(int argc, char** argv) -> int {
  auto options =
      cxxopts::Options("camel", "Camel is haplotype aware detection tool");

  /* clang-format off */
  options.add_options("serialization arguments")
    ("d,dst", "output destination folder",
      cxxopts::value<std::string>()->default_value("./camel_out"));
  options.add_options("correction arguments")
    ("o,n_overlaps", "number of overlaps per read to keep",
      cxxopts::value<std::size_t>()->default_value("64"))
    ("m,match", "score for matching bases",
      cxxopts::value<std::int8_t>()->default_value("3"))
    ("n,mismatch", "score for mismatching bases",
      cxxopts::value<std::int8_t>()->default_value("-5"))
    ("g,gap", "gap penalty (must be nagative)",
      cxxopts::value<std::int8_t>()->default_value("-4"))
    ("w,window_length", "targeted correction window len",
      cxxopts::value<std::uint32_t>()->default_value("500"))
    ("q,quality", "minimum read average quality",
      cxxopts::value<std::uint32_t>()->default_value("10"));
  options.add_options("utility arguments")
    ("t,threads", "number of threads avalable for execution",
            cxxopts::value<std::uint32_t>());
  options.add_options("input")
    ("overlaps", 
      "overlaps path", cxxopts::value<std::string>())
    ("reads", "input fastq reads", 
            cxxopts::value<std::vector<std::string>>());
  options.parse_positional({"overlaps", "reads"});
  /* clang-format on */

  auto const result = options.parse(argc, argv);

  auto const n_threads = result["threads"].as<std::uint32_t>();
  auto const dst_path = std::filesystem::path(result["dst"].as<std::string>());

  auto task_arena = tbb::task_arena(n_threads);

  if (std::filesystem::exists(dst_path)) {
    std::filesystem::remove_all(dst_path);
  }

  std::filesystem::create_directory(dst_path);

  auto read_paths = std::vector<std::filesystem::path>();
  auto overlaps_path =
      std::filesystem::path(result["overlaps"].as<std::string>());

  {
    auto paths_strs = result["reads"].as<std::vector<std::string>>();
    for (auto const& it : paths_strs) {
      auto it_path = std::filesystem::path(std::move(it));
      if (std::filesystem::is_regular_file(it_path)) {
        read_paths.emplace_back(std::move(it_path));
      } else if (std::filesystem::is_directory(it_path)) {
        for (auto dir_entry : std::filesystem::recursive_directory_iterator(
                 std::move(it_path))) {
          if (std::filesystem::is_regular_file(dir_entry)) {
            read_paths.emplace_back(std::move(dir_entry));
          }
        }
      }
    }
  }

  decltype(read_paths)(std::make_move_iterator(read_paths.begin()),
                       std::make_move_iterator(read_paths.end()))
      .swap(read_paths);

  auto timer = biosoup::Timer();

  task_arena.execute([&] {
    auto const correct_cfg = camel::CorrectConfig{
        .poa_cfg = camel::POAConfig{},
        .correct_window = result["window_length"].as<std::uint32_t>()};

    timer.Start();
    auto reads = camel::LoadSequences(read_paths);
    fmt::print(stderr, "[camel]({:12.3f}) loaded {} reads\n", timer.Stop(),
               reads.size());

    timer.Start();
    auto overlaps = camel::LoadOverlaps(overlaps_path, reads,
                                        result["n_overlaps"].as<std::size_t>());
    auto const n_ovlps = std::transform_reduce(
        overlaps.cbegin(), overlaps.cend(), 0ULL, std::plus<std::size_t>(),
        std::mem_fn(&std::vector<biosoup::Overlap>::size));

    fmt::print(stderr, "[camel]({:12.3f}) loaded {} overlaps\n", timer.Stop(),
               n_ovlps);

    timer.Start();
    task_arena.execute([&]() -> void {
      auto corrected_reads = camel::ErrorCorrect(correct_cfg, std::move(reads),
                                                 std::move(overlaps));
      camel::StoreSequences(corrected_reads, dst_path, 1U << 28U);
    });
    timer.Stop();
    fmt::print(stderr, "[camel]({:12.3f}) done\n", timer.elapsed_time());
  });

  return EXIT_SUCCESS;
}
