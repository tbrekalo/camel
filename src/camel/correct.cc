#include "camel/correct.h"

#include <array>
#include <functional>
#include <numeric>
#include <random>
#include <variant>

#include "biosoup/timer.hpp"
#include "detail/coverage.h"
#include "detail/window.h"
#include "fmt/core.h"
#include "tbb/parallel_for.h"

namespace camel {

auto ErrorCorrect(CorrectConfig const correct_cfg,
                  std::vector<std::unique_ptr<biosoup::NucleicAcid>> reads,
                  std::vector<std::vector<biosoup::Overlap>> overlaps)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
  auto corrected_targets = std::vector<std::unique_ptr<biosoup::NucleicAcid>>();
  auto dst = std::vector<std::unique_ptr<biosoup::NucleicAcid>>();
  auto function_timer = biosoup::Timer();

  function_timer.Start();
  auto coverage_estimate = detail::EstimateCoverage(reads, overlaps);
  fmt::print(stderr, "[camel::ErrorCorrect]({:12.3f}) coverage estimate {}\n",
             function_timer.Stop(), coverage_estimate);

  function_timer.Start();
  auto target_ids = std::vector<std::size_t>();
  for (auto i = 0U; i < overlaps.size(); ++i) {
    if (!overlaps[i].empty()) {
      target_ids.push_back(i);
    }
  }

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
      fmt::print(stderr,
                 "\r[camel::ErrorCorrect]({:12.3f}) aligned {:3.3f}% | polished "
                 "{:3.3f}%",
                 function_timer.Lap(), to_percent(n_aligned),
                 to_percent(n_polished));
    }
  };

  dst.resize(n_targets);

  function_timer.Start();
  tbb::parallel_for(
      std::size_t(0), target_ids.size(), [&](std::size_t target_index) -> void {
        auto const read_id = target_ids[target_index];
        auto alignments =
            std::vector<EdlibAlignResult>(overlaps[read_id].size());
        tbb::parallel_for(std::size_t(0), overlaps[read_id].size(),
                          [&, read_id](std::size_t overlap_id) -> void {
                            alignments[overlap_id] = detail::OverlapToALignment(
                                reads, overlaps[read_id][overlap_id]);
                          });
        ++n_aligned;
        report_state();

        auto windows = detail::CreateWindowsFromAlignments(
            reads, overlaps[read_id], alignments, correct_cfg.correct_window,
            coverage_estimate);

        auto backbone = reads[read_id]->InflateData();
        auto window_consensus = std::vector<std::string>(windows.size());

        tbb::parallel_for(
            std::size_t(0), windows.size(), [&](std::size_t window_id) -> void {
              auto interval = windows[window_id].interval;
              window_consensus[window_id] = detail::WindowConsensus(
                  std::string_view(std::next(backbone.cbegin(), interval.first),
                                   std::next(backbone.cbegin(), interval.last)),
                  detail::ReferenceWindowView{
                      .interval = interval,
                      .aligned_segments = windows[window_id].aligned_segments},
                  correct_cfg.poa_cfg);
            });

        {
          auto consensus = std::string();
          consensus.reserve(reads[read_id]->inflated_len * 1.2);
          auto prev = 0U;
          for (auto win_idx = 0U; win_idx < windows.size(); ++win_idx) {
            auto const& interval = windows[win_idx].interval;

            consensus.insert(consensus.end(), 
                             std::next(backbone.begin(), prev),
                             std::next(backbone.begin(), interval.first));
            consensus += window_consensus[win_idx];
            prev = interval.last;
          }

          consensus.insert(consensus.end(), std::next(backbone.begin(), prev),
                           backbone.end());

          dst[target_index] = std::make_unique<biosoup::NucleicAcid>(
              fmt::format("polished_{:09d}", n_polished), consensus);
        }

        ++n_polished;
        report_state();

        overlaps[read_id].clear();
      });

  detail::ReleaseAlignmentEngines();

  report_state();
  fmt::print(stderr, "\n");

  fmt::print(stderr, "[camel::ErrorCorrect]({:12.3f})\n",
             function_timer.elapsed_time());
  return dst;
}

}  // namespace camel
