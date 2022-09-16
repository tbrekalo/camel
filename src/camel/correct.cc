#include "camel/correct.h"

#include <array>
#include <functional>
#include <numeric>
#include <random>
#include <variant>

#include "biosoup/timer.hpp"
#include "detail/coverage.h"
#include "detail/overload.h"
#include "detail/task_queue.h"
#include "fmt/core.h"
#include "tsl/robin_map.h"

namespace camel {

namespace detail {

static auto constexpr kUpdateInterval = 1. / 24.;

struct ConsensuState {
  std::string backbone;
  std::vector<ReferenceWindow> ref_windows;
  std::vector<std::string> window_consensus;

  std::size_t n_awaiting_tasks;
  tsl::robin_map<TaskIdType, std::uint32_t> task_mapping;
};

}  // namespace detail

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
  auto const n_targets = std::transform_reduce(
      overlaps.cbegin(), overlaps.cend(), 0UL, std::plus<size_t>(),
      [](std::vector<biosoup::Overlap> const& ovlp_vec) -> size_t {
        return !ovlp_vec.empty();
      });

  fmt::print(stderr, "[camel::ErrorCorrect]({:12.3f}) n target reads {}\n",
             function_timer.Stop(), n_targets);

  auto n_aligned = 0U;
  auto n_windowed = 0U;
  auto n_polished = 0U;
  auto n_fails = 0U;

  auto const report_status = [&]() -> void {
    auto const to_percent = [n_targets](double n) -> double {
      return 100. * n / n_targets;
    };
    fmt::print(stderr,
               "\r[camel::ErrorCorrect]({:12.3f}) aligned : {:3.2f}% | "
               "converted to windows : {:3.2f}% | "
               "polished : {:3.2f}% | failed polishes : {:3.2f}%",
               function_timer.Lap(), to_percent(n_aligned),
               to_percent(n_windowed), to_percent(n_polished),
               to_percent(n_fails));
  };

  auto task_queue = detail::TaskQueue();

  auto alignment_registry =
      tsl::robin_map<detail::TaskIdType,
                     std::pair<std::uint32_t, std::uint32_t>>(
          reads.size() * coverage_estimate);
  auto alignment_counts = std::vector<std::uint32_t>(reads.size());

  auto windowing_registry = tsl::robin_map<detail::TaskIdType, std::uint32_t>(
      reads.size() * coverage_estimate);
  auto consensus_registry = tsl::robin_map<detail::TaskIdType, std::uint32_t>(
      reads.size() * coverage_estimate);

  auto alignments = std::vector<std::vector<EdlibAlignResult>>(reads.size());
  auto consensus_states = std::vector<detail::ConsensuState>(reads.size());

  for (auto i = 0U; i < reads.size(); ++i) {
    alignments[i].resize(overlaps[i].size());
    alignment_counts[i] = overlaps[i].size();
  }

  for (auto read_idx = 0U; read_idx < reads.size(); ++read_idx) {
    if (!overlaps[read_idx].empty()) {
      for (auto ovlp_idx = 0U; ovlp_idx < overlaps[read_idx].size();
           ++ovlp_idx) {
        alignment_registry[task_queue.Push(
            detail::AlignmentArgPack(reads, overlaps[read_idx][ovlp_idx]))] =
            std::pair(read_idx, ovlp_idx);
      }
    } else {
      dst.push_back(std::make_unique<biosoup::NucleicAcid>(*reads[read_idx]));
    }
  }

  function_timer.Start();
  {
    auto update_timer = biosoup::Timer();

    update_timer.Start();
    while (n_polished + n_fails < n_targets) {
      if (auto opt_result = task_queue.TryPop(); opt_result) {
        std::visit(
            detail::overload{
                [&, task_id = opt_result->task_id](
                    detail::AlignmentResult alignment) -> void {
                  auto const [read_id, ovlp_id] = alignment_registry[task_id];

                  alignments[read_id][ovlp_id] = alignment;
                  if (--alignment_counts[read_id] == 0U) {
                    windowing_registry[task_queue.Push(
                        detail::WindowConstructArgPack(
                            reads, overlaps[read_id], alignments[read_id],
                            correct_cfg.correct_window, coverage_estimate))] =
                        read_id;
                    if (++n_aligned == n_targets) {
                      std::remove_cvref_t<decltype(alignment_registry)>().swap(
                          alignment_registry);
                      std::remove_cvref_t<decltype(alignment_counts)>().swap(
                          alignment_counts);
                    }
                  }
                },
                [&, task_id = opt_result->task_id](
                    detail::WindowConstructResult ref_windows) -> void {
                  auto const read_id = windowing_registry[task_id];
                  auto& consensus_state = consensus_states[read_id];
                  auto const n_windows = ref_windows.size();

                  if (n_windows != 0) {
                    {
                      consensus_state = detail::ConsensuState{
                          .backbone = reads[read_id]->InflateData(),
                          .ref_windows = std::move(ref_windows),
                          .window_consensus =
                              std::vector<std::string>(n_windows),
                          .n_awaiting_tasks = n_windows,
                          .task_mapping =
                              tsl::robin_map<detail::TaskIdType, std::uint32_t>(
                                  n_windows)};
                    }

                    for (auto win_id = 0U;
                         win_id < consensus_state.ref_windows.size();
                         ++win_id) {
                      auto& ref_win = consensus_state.ref_windows[win_id];
                      auto const win_task_id =
                          task_queue.Push(detail::WindowConsensusArgPack(
                              std::string_view(
                                  std::next(consensus_state.backbone.begin(),
                                            ref_win.interval.first),
                                  std::next(consensus_state.backbone.begin(),
                                            ref_win.interval.last)),
                              detail::ReferenceWindowView{
                                  .interval = ref_win.interval,
                                  .aligned_segments = ref_win.aligned_segments},
                              correct_cfg.poa_cfg));

                      consensus_state.task_mapping[win_task_id] = win_id;
                      consensus_registry[win_task_id] = read_id;
                    }

                  } else {
                    dst.push_back(std::make_unique<biosoup::NucleicAcid>(
                        *reads[read_id]));
                    ++n_fails;
                  }

                  ++n_windowed;

                  std::vector<biosoup::Overlap>{}.swap(overlaps[read_id]);
                  std::vector<EdlibAlignResult>{}.swap(alignments[read_id]);

                  if (n_windowed + n_fails == n_targets) {
                    std::remove_cvref_t<decltype(windowing_registry)>().swap(
                        windowing_registry);
                  }
                },
                [&, task_id = opt_result->task_id](
                    detail::WindowConsensusResult win_consensus) -> void {
                  auto const read_id = consensus_registry[task_id];
                  auto& consensus_state = consensus_states[read_id];

                  consensus_state
                      .window_consensus[consensus_state.task_mapping[task_id]] =
                      std::move(win_consensus);

                  if (--consensus_state.n_awaiting_tasks == 0U) {
                    auto consensus = std::string();
                    consensus.reserve(reads[read_id]->inflated_len * 1.2);
                    auto prev = 0U;
                    for (auto win_idx = 0U;
                         win_idx < consensus_state.window_consensus.size();
                         ++win_idx) {
                      auto const& interval =
                          consensus_state.ref_windows[win_idx].interval;

                      consensus.insert(
                          consensus.end(),
                          std::next(consensus_state.backbone.begin(), prev),
                          std::next(consensus_state.backbone.begin(),
                                    interval.first));
                      consensus += consensus_state.window_consensus[win_idx];
                      prev = interval.last;
                    }

                    consensus.insert(
                        consensus.end(),
                        std::next(consensus_state.backbone.begin(), prev),
                        consensus_state.backbone.end());

                    dst.emplace_back(std::make_unique<biosoup::NucleicAcid>(
                        fmt::format("polished_{:09d}", n_polished), consensus));

                    consensus_state = detail::ConsensuState{};
                    ++n_polished;
                  }
                }},
            std::move(opt_result->value));
      }

      if (update_timer.Lap() > detail::kUpdateInterval) {
        report_status();
        update_timer.Start();
      }
    }

    report_status();
    fmt::print(stderr, "\n");
    fmt::print(stderr, "[camel::ReleaseAlignmentEngines] released {} engines\n",
               detail::ReleaseAlignmentEngines());
  }

  fmt::print(stderr, "[camel::ErrorCorrect]({:12.3f})\n",
             function_timer.elapsed_time());
  return dst;
}

}  // namespace camel
