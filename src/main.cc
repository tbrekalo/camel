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
    ("o,out", "output destination folder for reads.fa",
      cxxopts::value<std::string>()->default_value("./"));
  options.add_options("correction arguments")
    ("n,n-overlaps", "number of overlaps per read to keep",
      cxxopts::value<std::size_t>()->default_value("64"))
    ("e,error-threshold", "maximum allowed overlap error threshold",
      cxxopts::value<double>()->default_value("0.3"))
    ("m,match", "score for matching bases",
      cxxopts::value<std::int8_t>()->default_value("3"))
    ("x,mismatch", "score for mismatching bases",
      cxxopts::value<std::int8_t>()->default_value("-5"))
    ("g,gap", "gap penalty (must be nagative)",
      cxxopts::value<std::int8_t>()->default_value("-4"))
    ("w,window-length", "targeted correction window len",
      cxxopts::value<std::uint32_t>()->default_value("200"))
    ("q,quality", "minimum read average quality",
      cxxopts::value<std::uint32_t>()->default_value("10"));
  options.add_options("utility arguments")
    ("t,threads", "number of threads avalable for execution",
            cxxopts::value<std::uint32_t>()->default_value("1"));
  options.add_options("input")
    ("reads", "input fastq reads", 
            cxxopts::value<std::string>())
    ("overlaps", "overlaps path", 
            cxxopts::value<std::string>());
  options.add_options()
    ("h,help", "print help");
  options.parse_positional({"reads", "overlaps"});
  /* clang-format on */

  try {
    auto const result = options.parse(argc, argv);
    if (result.count("help")) {
      fmt::print(stderr, "{}", options.help());
      return EXIT_SUCCESS;
    }

    auto const n_threads = result["threads"].as<std::uint32_t>();
    auto const out_path =
        std::filesystem::path(result["out"].as<std::string>());

    auto task_arena = tbb::task_arena(n_threads);
    auto reads_path = std::filesystem::path(result["reads"].as<std::string>());
    auto overlaps_path =
        std::filesystem::path(result["overlaps"].as<std::string>());

    auto timer = biosoup::Timer();

    task_arena.execute([&] {
      auto const correct_cfg = camel::CorrectConfig{
          .poa_cfg = camel::POAConfig{},
          .window_length = result["window-length"].as<std::uint32_t>()};

      timer.Start();
      auto reads = camel::LoadSequences(reads_path);
      fmt::print(stderr, "[camel]({:12.3f}) loaded {} reads\n", timer.Stop(),
                 reads.size());

      timer.Start();
      auto overlaps = camel::LoadOverlaps(
          overlaps_path, reads, result["error-threshold"].as<double>(),
          result["n-overlaps"].as<std::size_t>());
      auto const n_ovlps = std::transform_reduce(
          overlaps.cbegin(), overlaps.cend(), 0ULL, std::plus<std::size_t>(),
          std::mem_fn(&std::vector<biosoup::Overlap>::size));

      fmt::print(stderr, "[camel]({:12.3f}) loaded {} overlaps\n", timer.Stop(),
                 n_ovlps);

      timer.Start();
      task_arena.execute([&]() -> void {
        auto corrected_reads = camel::ErrorCorrect(
            correct_cfg, std::move(reads), std::move(overlaps));
        camel::StoreSequences(corrected_reads, out_path);
      });
      timer.Stop();
      fmt::print(stderr, "[camel]({:12.3f}) done\n", timer.elapsed_time());
    });
  } catch (std::exception const& e) {
    fmt::print(stderr, "{}", e.what());
  }

  return EXIT_SUCCESS;
}
