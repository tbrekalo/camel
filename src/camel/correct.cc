#include "camel/correct.h"

#include <array>
#include <functional>
#include <iterator>
#include <numeric>
#include <random>
#include <variant>

#include "biosoup/timer.hpp"
#include "detail/coverage.h"
#include "detail/window.h"
#include "fmt/core.h"
#include "tbb/parallel_for.h"

namespace camel {

namespace detail {

static auto GetTargetIds(
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> reads,
    std::span<std::vector<biosoup::Overlap> const> overlaps)
    -> std::vector<std::size_t> {
  auto dst = std::vector<std::size_t>();
  for (auto i = 0U; i < overlaps.size(); ++i) {
    if (!overlaps[i].empty()) {
      dst.push_back(i);
    }
  }

  return dst;
}

static auto AlignOverlaps(
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> reads,
    std::span<biosoup::Overlap> overlaps) -> void {
  tbb::parallel_for(std::size_t(0), overlaps.size(),
                    [reads, overlaps](std::size_t overlap_id) -> void {
                      overlaps[overlap_id] =
                          AlignedOverlap(reads, overlaps[overlap_id]);
                    });
}

static auto CreateConsensus(POAConfig poa_cfg, NucleicView read,
                            std::span<ReferenceWindow const> windows) {
  auto backbone_data = read.InflateData();
  auto backbone_quality = read.InflateQuality();

  auto const fetch_quality = [&](std::uint32_t first, std::uint32_t last) {
    if (!backbone_quality.empty()) {
      std::string_view(std::next(backbone_quality.cbegin(), first),
                       std::next(backbone_quality.cbegin(), last));
    }

    return std::string_view{};
  };

  auto consensuses = std::vector<ConsensusResult>(windows.size());

  tbb::parallel_for(
      std::size_t(0), windows.size(),
      [poa_cfg, windows, &backbone_data, &fetch_quality,
       &consensuses](std::size_t window_id) -> void {
        auto interval = windows[window_id].interval;
        consensuses[window_id] = WindowConsensus(
            std::string_view(std::next(backbone_data.cbegin(), interval.first),
                             std::next(backbone_data.cbegin(), interval.last)),
            fetch_quality(interval.first, interval.last),
            ReferenceWindowView{
                .interval = interval,
                .aligned_segments = windows[window_id].aligned_segments},
            poa_cfg);
      });

  auto consensus = std::string();
  consensus.reserve(read.InflatedLenght() * 1.2);
  for (auto win_idx = 0U; win_idx < windows.size(); ++win_idx) {
    consensus += consensuses[win_idx].bases;
  }

  auto const polished_ratio =
      std::transform_reduce(consensuses.cbegin(), consensuses.cend(), 0.0,
                            std::plus<double>(),
                            [](detail::ConsensusResult const& cr) {
                              return static_cast<double>(cr.is_corrected);
                            }) /
      static_cast<double>(consensuses.size());

  /* clang-format off */
          std::string consensus_name =
              fmt::format("{} LN:i:{} RC:i:{} XC:f:{:5.2f}",
                  read.Name(),
                          consensus.size(), windows.size(), polished_ratio);
  /* clang-format on */

  return std::make_unique<biosoup::NucleicAcid>(consensus_name, consensus);
}

}  // namespace detail

auto ErrorCorrect(CorrectConfig const correct_cfg,
                  std::vector<std::unique_ptr<biosoup::NucleicAcid>> reads,
                  std::vector<std::vector<biosoup::Overlap>> overlaps)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
  auto dst = std::vector<std::unique_ptr<biosoup::NucleicAcid>>();
  auto function_timer = biosoup::Timer();

  function_timer.Start();
  auto coverage_estimate = detail::EstimateCoverage(reads, overlaps);
  fmt::print(stderr, "[camel::ErrorCorrect]({:12.3f}) coverage estimate {}\n",
             function_timer.Stop(), coverage_estimate);

  function_timer.Start();

  auto const target_ids = detail::GetTargetIds(reads, overlaps);
  auto const n_targets = target_ids.size();
  fmt::print(stderr, "[camel::ErrorCorrect]({:12.3f}) n target reads {}\n",
             function_timer.Stop(), n_targets);

  auto n_aligned = std::atomic_size_t(0);
  auto n_polished = std::atomic_size_t(0);
  auto report_ticket = std::atomic_size_t(0);
  auto const report_state = [&function_timer, n_targets, &n_aligned,
                             &n_polished, &report_ticket]() -> void {
    auto const to_percent = [n_targets](double val) -> double {
      return 100. * (val / n_targets);
    };

    if (auto ticket = ++report_ticket; ticket == report_ticket) {
      fmt::print(
          stderr,
          "\r[camel::ErrorCorrect]({:12.3f}) aligned {:3.3f}% | polished "
          "{:3.3f}%",
          function_timer.Lap(), to_percent(n_aligned), to_percent(n_polished));
    }
  };

  dst.resize(n_targets);

  function_timer.Start();
  tbb::parallel_for(
      std::size_t(0), target_ids.size(), [&](std::size_t target_idx) -> void {
        auto const read_id = target_ids[target_idx];
        detail::AlignOverlaps(reads, overlaps[read_id]);
        ++n_aligned;
        report_state();

        auto windows = detail::CreateWindowsFromAlignments(
            reads, overlaps[read_id], correct_cfg.window_cfg,
            coverage_estimate);

        dst[target_idx] = detail::CreateConsensus(
            correct_cfg.poa_cfg,
            detail::NucleicView(reads[read_id].get(), false), windows);

        ++n_polished;
        report_state();

        overlaps[read_id].clear();
      });

  detail::ReleaseAlignmentEngines();

  report_state();
  fmt::print(stderr, "\n");

  function_timer.Stop();
  fmt::print(stderr, "[camel::ErrorCorrect]({:12.3f})\n",
             function_timer.elapsed_time());
  return dst;
}

}  // namespace camel
