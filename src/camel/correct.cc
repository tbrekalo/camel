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

  auto n_aligned = 0;
  auto n_windowed = 0;
  auto n_polished = 0;

  auto const report_status = [&]() -> void {
    auto const to_percent = [n_targets](double n) -> double {
      return 100. * n / n_targets;
    };
    fmt::print(stderr,
               "\r[camel::ErrorCorrect]({:12.3f}) aligned : {:3.2f}% | "
               "n_windowed : {:3.2f}% | "
               "n_polished : {:3.2f}%",
               function_timer.Lap(), to_percent(n_aligned),
               to_percent(n_windowed), to_percent(n_polished));
  };

  auto task_queue = detail::TaskQueue();

  auto alignment_registry =
      tsl::robin_map<detail::TaskIdType,
                     std::pair<std::uint32_t, std::uint32_t>>(
          reads.size() * coverage_estimate);
  auto alignment_counts = std::vector<std::uint32_t>(reads.size());
  auto windowing_registry = tsl::robin_map<detail::TaskIdType, std::uint32_t>(
      reads.size() * coverage_estimate);

  auto alignments = std::vector<std::vector<EdlibAlignResult>>(reads.size());
  auto ref_windows =
      std::vector<std::vector<detail::ReferenceWindow>>(reads.size());

  for (auto i = 0U; i < reads.size(); ++i) {
    alignments[i].resize(overlaps[i].size());
    alignment_counts[i] = overlaps[i].size();
  }

  for (auto read_idx = 0U; read_idx < reads.size(); ++read_idx) {
    for (auto ovlp_idx = 0U; ovlp_idx < overlaps[read_idx].size(); ++ovlp_idx) {
      alignment_registry[task_queue.Push(
          detail::AlignmentArgPack(reads, overlaps[read_idx][ovlp_idx]))] =
          std::pair(read_idx, ovlp_idx);
    }
  }

  function_timer.Start();
  {
    auto update_timer = biosoup::Timer();

    update_timer.Start();
    while (n_windowed < n_targets) {
      if (auto opt_result = task_queue.TryPop(); opt_result) {
        std::visit(
            detail::overload{
                [&, task_id = opt_result->task_id](
                    detail::AlignmentResult alignment) -> void {
                  auto const [read_id, ovlp_id] = alignment_registry[task_id];

                  alignments[read_id][ovlp_id] = alignment;
                  if (--alignment_counts[read_id] == 0U) {
                    windowing_registry[task_queue.Push(detail::WindowArgPack(
                        reads, overlaps[read_id], alignments[read_id],
                        correct_cfg.correct_window, coverage_estimate))] =
                        read_id;

                    ++n_aligned;
                  }
                },
                [&, task_id = opt_result->task_id](
                    detail::WindowResult ref_windows) -> void {
                  auto const read_id = windowing_registry[task_id];
                  ++n_windowed;
                },
                [](auto&&) -> void {}},
            std::move(opt_result->value));
      }

      if (update_timer.Lap() > detail::kUpdateInterval) {
        report_status();
        update_timer.Start();
      }
    }

    report_status();
    fmt::print(stderr, "\n");

    // for (auto const& it : alignments) {
    //   for (auto const& jt : it) {
    //     edlibFreeAlignResult(jt);
    //   }
    // }
  }

  // auto windowing_registry =

  // auto target_ids = std::vector<std::uint32_t>();
  // for (auto read_id = 0U; read_id < overlaps.size(); ++read_id) {
  //   if (!overlaps[read_id].empty()) {
  //     target_ids.push_back(read_id);
  //   }
  // }

  // auto const kNTargets = target_ids.size();
  // auto const kNSupporting = overlaps.size() - kNTargets;

  // auto const kNOverlaps = std::transform_reduce(
  //     overlaps.cbegin(), overlaps.cend(), 0UL, std::plus<std::size_t>(),
  //     [](std::vector<biosoup::Overlap> const& ovlps) -> std::size_t {
  //       return ovlps.size();
  //     });

  // fmt::print(
  //     stderr,
  //     "[camel::ErrorCorrect]({:12.3f}) (kNTargets, kNSupporting) = ({},
  //     {})\n", timer.Stop(), kNTargets, kNSupporting);

  // timer.Start();

  // {
  //   timer.Start();
  //   auto align_futures =
  //       std::vector<std::vector<std::future<EdlibAlignResult>>>(kNTargets);

  //   for (auto i = 0U; i < kNTargets; ++i) {
  //     auto const target_id = target_ids[i];
  //     align_futures[i].resize(overlaps[target_id].size());
  //     for (auto j = 0U; j < overlaps[target_id].size(); ++j) {
  //       align_futures[i][j] = state.thread_pool->Submit(
  //           [&src_reads, &ovlp = overlaps[target_id][j]]() ->
  //           EdlibAlignResult {
  //             auto const [lhs_str, rhs_str] =
  //                 detail::ExtractSubstrings(src_reads, ovlp);
  //             return detail::AlignStrings(lhs_str, rhs_str);
  //           });
  //     }
  //   }

  //   auto window_futures =
  //       std::vector<std::future<std::vector<detail::ReferenceWindow>>>(
  //           align_futures.size());

  //   auto n_transformed = 0U;
  //   for (auto i = 0U; i < kNTargets; ++i) {
  //     auto alignments =
  //     std::vector<EdlibAlignResult>(align_futures[i].size());
  //     std::transform(std::make_move_iterator(align_futures[i].begin()),
  //                    std::make_move_iterator(align_futures[i].end()),
  //                    alignments.begin(),
  //                    std::mem_fn(&std::future<EdlibAlignResult>::get));

  //     n_transformed += alignments.size();
  //     window_futures[i] = state.thread_pool->Submit(
  //         [&src_reads, &overlaps, kCovgEstimate,
  //          alignments = std::move(alignments)](
  //             std::uint32_t read_id) ->
  //             std::vector<detail::ReferenceWindow>
  //             {
  //           auto windows = detail::CreateWindowsFromAlignments(
  //               src_reads, std::move(overlaps[read_id]),
  //               std::move(alignments), kCovgEstimate);

  //           return windows;
  //         },
  //         target_ids[i]);

  //     fmt::print(stderr,
  //                "\r[camel::ErrorCorrect({:12.3f})] transformed {} / {} "
  //                "alignment futures to window tasks",
  //                timer.Lap(), n_transformed, kNOverlaps);
  //   }

  //   fmt::print(stderr,
  //              "\r[camel::ErrorCorrect]({:12.3f}) transformed {} / {} "
  //              "alignment futures "
  //              "to window tasks\n",
  //              timer.Lap(), n_transformed, kNOverlaps);
  //   decltype(align_futures){}.swap(align_futures);

  //   auto ref_windows =
  //   std::vector<std::vector<detail::ReferenceWindow>>();
  //   ref_windows.reserve(window_futures.size());

  //   for (auto i = 0U; i < window_futures.size(); ++i) {
  //     ref_windows.push_back(window_futures[i].get());

  //     if ((i & 127) == 0U) {
  //       fmt::print(stderr,
  //                  "\r[camel::ErrorCorrect]({:12.3f}) collected window
  //                  futures "
  //                  "{} / {}",
  //                  timer.Lap(), i + 1, window_futures.size());
  //     }
  //   }

  //   fmt::print(stderr,
  //              "\r[camel::ErrorCorrect]({:12.3f}) collected window
  //              futures "
  //              "{} / {}\n",
  //              timer.Stop(), window_futures.size(),
  //              window_futures.size());
  //   decltype(window_futures)().swap(window_futures);

  //   auto alignment_engines =
  //       tsl::robin_map<std::thread::id,
  //                      std::unique_ptr<spoa::AlignmentEngine>>();

  //   for (auto const [thread_id, _] : state.thread_pool->thread_map()) {
  //     alignment_engines[thread_id] = spoa::AlignmentEngine::Create(
  //         spoa::AlignmentType::kNW, correct_cfg.poa_cfg.match,
  //         correct_cfg.poa_cfg.mismatch, correct_cfg.poa_cfg.gap);
  //   }

  //   timer.Start();
  //   dst.reserve(window_futures.size());
  //   auto spoa_futures = std::vector<std::future<void>>();

  //   for (auto read_idx = 0U; read_idx < ref_windows.size(); ++read_idx) {
  //     auto graphs =
  //     std::vector<spoa::Graph>(ref_windows[read_idx].size()); auto
  //     backbone = src_reads[target_ids[read_idx]]->InflateData();

  //     for (auto win_idx = 0U; win_idx < ref_windows[read_idx].size();
  //          ++win_idx) {
  //       graphs[win_idx].AddAlignment(
  //           spoa::Alignment(),
  //           backbone.substr(ref_windows[read_idx][win_idx].interval.first,
  //                           detail::IntervalLength(
  //                               ref_windows[read_idx][win_idx].interval)));
  //     }

  //     for (auto win_idx = 0U; win_idx < ref_windows[read_idx].size();
  //          ++win_idx) {
  //       spoa_futures.emplace_back(state.thread_pool->Submit(
  //           [&alignment_engines,
  //            active_window = std::move(ref_windows[read_idx][win_idx]),
  //            &graphs, win_idx]() -> void {
  //             auto const& alignment_engine =
  //                 alignment_engines[std::this_thread::get_id()];

  //             auto const& [win_ref_intv, aligned_segments] =
  //             active_window; auto const window_length =
  //             IntervalLength(win_ref_intv);

  //             for (auto const& [alignment_interval, bases] :
  //             aligned_segments) {
  //               auto const alignment_len =
  //               IntervalLength(alignment_interval); if (alignment_len <
  //                   window_length * detail::kSmallWindowPercent) {
  //                 continue;
  //               }

  //               auto const legal_start =
  //                   window_length * detail::kAllowedFuzzPercent;
  //               auto const legal_end = window_length - legal_start;

  //               auto alignment = spoa::Alignment();
  //               if (alignment_interval.first <= legal_start &&
  //                   legal_end <= alignment_interval.last) {
  //                 alignment = alignment_engine->Align(bases,
  //                 graphs[win_idx]);
  //               } else {
  //                 auto mapping = std::vector<spoa::Graph::Node const*>();
  //                 auto subgraph = graphs[win_idx].Subgraph(
  //                     alignment_interval.first, alignment_interval.last -
  //                     1, &mapping);

  //                 alignment =
  //                     alignment_engines[std::this_thread::get_id()]->Align(
  //                         bases, subgraph);
  //                 subgraph.UpdateAlignment(mapping, &alignment);
  //               }

  //               graphs[win_idx].AddAlignment(alignment, bases);
  //             }
  //           }));
  //     }

  //     auto consensus = std::string();
  //     consensus.reserve(src_reads[target_ids[read_idx]]->inflated_len
  //     * 1.2);

  //     for (auto i = 0; i < spoa_futures.size(); ++i) {
  //       spoa_futures[i].wait();
  //     }

  //     {
  //       auto prev = 0U;
  //       for (auto win_idx = 0U; win_idx < spoa_futures.size(); ++win_idx)
  //       {
  //         spoa_futures[win_idx].wait();
  //         auto const& interval = ref_windows[read_idx][win_idx].interval;

  //         consensus.insert(consensus.end(), std::next(backbone.begin(),
  //         prev),
  //                          std::next(backbone.begin(), interval.first));
  //         consensus += graphs[win_idx].GenerateConsensus();
  //         prev = interval.last;
  //       }

  //       consensus.insert(consensus.end(), std::next(backbone.begin(),
  //       prev),
  //                        backbone.end());
  //     }

  //     dst.push_back(std::make_unique<biosoup::NucleicAcid>(
  //         src_reads[target_ids[read_idx]]->name, consensus));
  //     fmt::print(
  //         stderr,
  //         "\r[camel::ErrorCorrect]({:12.3f}) generated consensus for {} /
  //         "
  //         "{} reads",
  //         timer.Lap(), read_idx + 1, ref_windows.size());
  //     spoa_futures.clear();
  //   }
  //   fmt::print(stderr,
  //              "\r[camel::ErrorCorrect]({:12.3f}) generated consensus for
  //              {} / "
  //              "{} reads\n",
  //              timer.Stop(), ref_windows.size(), ref_windows.size());
  // }

  fmt::print(stderr, "[camel::ErrorCorrect]({:12.3f})\n",
             function_timer.elapsed_time());
  return dst;
}

}  // namespace camel
