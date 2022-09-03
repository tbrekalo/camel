#ifndef CAMEL_DETAIL_TASK_QUEUE_H_
#define CAMEL_DETAIL_TASK_QUEUE_H_

#include <optional>
#include <tuple>
#include <variant>

#include "window.h"

namespace camel::detail {

using TaskIdType = std::uint64_t;

using AlignmentResult = EdlibAlignResult;
using WindowResult = std::vector<ReferenceWindow>;
using ConsensusResult = std::unique_ptr<biosoup::NucleicAcid>;

using ResultVariant =
    std::variant<AlignmentResult, WindowResult, ConsensusResult>;

struct TaskResult {
  TaskIdType task_id;
  ResultVariant value;
};

using AlignmentArgPack =
    std::tuple<nonstd::span<std::unique_ptr<biosoup::NucleicAcid>>,
               biosoup::Overlap>;
using WindowArgPack =
    std::tuple<nonstd::span<std::unique_ptr<biosoup::NucleicAcid>>,
               nonstd::span<biosoup::Overlap>, nonstd::span<EdlibAlignResult>,
               std::uint32_t>;

using ArgsPack = std::variant<AlignmentArgPack, WindowArgPack>;

class TaskQueue {
 public:
  explicit TaskQueue();

  [[nodiscard]] auto Push(ArgsPack args) -> TaskIdType;
  [[nodiscard]] auto TryPop() -> std::optional<TaskResult>;

 private:
  struct Impl;
  std::shared_ptr<Impl> pimpl_;
};

}  // namespace camel::detail

#endif /* CAMEL_DETAIL_TASK_QUEUE_H_ */
