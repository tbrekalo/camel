#include <filesystem>
#include <iostream>

#include "biosoup/timer.hpp"
#include "camel/correct.h"
#include "camel/io.h"
#include "camel/version.h"
#include "cxxopts.hpp"
#include "fmt/core.h"
#include "tbb/task_arena.h"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects(0);

auto main(int argc, char** argv) -> int {
  auto options = cxxopts::Options("camel", "Camel is a read correction tool");

  /* clang-format off */
  options.add_options("correction arguments")
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
    ("q,quality_threshold", "minimum read average quality",
      cxxopts::value<double>()->default_value("10.0"));
  options.add_options("utility arguments")
    ("t,threads", "number of threads avalable for execution",
            cxxopts::value<std::uint32_t>()->default_value("1"));
  options.add_options("input")
    ("reads", "input fastq reads", 
            cxxopts::value<std::string>())
    ("overlaps", "overlaps path", 
            cxxopts::value<std::string>());
  options.positional_help("<reads_path> <overlaps_path>");
  options.add_options("info")
    ("h,help", "print help")
    ("v,version", "print version and exit");
  /* clang-format on */

  try {
    options.show_positional_help();
    options.parse_positional({"reads", "overlaps"});
    auto const result = options.parse(argc, argv);
    if (result.count("help")) {
      fmt::print(stderr, "{}", options.help());
      return EXIT_SUCCESS;
    }

    /* clang-format off */
    if (result.count("version")) {
      fmt::print(stderr, "{}.{}.{}", 
          camel_VERSION_MAJOR, 
          camel_VERSION_MINOR,
          camel_VERSION_PATCH);
      return EXIT_SUCCESS;
    }
    /* clang-format on */

    auto const n_threads = result["threads"].as<std::uint32_t>();

    auto task_arena = tbb::task_arena(n_threads);
    auto reads_path = std::filesystem::path(result["reads"].as<std::string>());
    auto overlaps_path =
        std::filesystem::path(result["overlaps"].as<std::string>());

    auto timer = biosoup::Timer();

    task_arena.execute([&] {
      auto const correct_cfg = camel::CorrectConfig{
          .poa_cfg = camel::POAConfig{},
          .window_cfg = {
              .quality_threshold = result["quality_threshold"].as<double>(),
              .window_length = result["window-length"].as<std::uint32_t>()}};

      timer.Start();
      auto reads = camel::LoadSequences(reads_path);
      fmt::print(stderr, "[camel]({:12.3f}) loaded {} reads\n", timer.Stop(),
                 reads.size());

      timer.Start();
      auto overlaps = camel::LoadOverlaps(
          overlaps_path, reads, result["error-threshold"].as<double>());
      auto const n_ovlps = std::transform_reduce(
          overlaps.cbegin(), overlaps.cend(), 0ULL, std::plus<std::size_t>(),
          std::mem_fn(&std::vector<biosoup::Overlap>::size));

      fmt::print(stderr, "[camel]({:12.3f}) loaded {} overlaps\n", timer.Stop(),
                 n_ovlps);

      timer.Start();
      task_arena.execute([&]() -> void {
        auto corrected_reads = camel::ErrorCorrect(
            correct_cfg, std::move(reads), std::move(overlaps));
        camel::SerializeSequences(corrected_reads, std::cout);
      });
      timer.Stop();
      fmt::print(stderr, "[camel]({:12.3f}) done\n", timer.elapsed_time());
    });
  } catch (std::exception const& e) {
    fmt::print(stderr, "{}", e.what());
  }

  return EXIT_SUCCESS;
}
