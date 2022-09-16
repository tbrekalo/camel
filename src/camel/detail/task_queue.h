#ifndef CAMEL_DETAIL_TASK_QUEUE_H_
#define CAMEL_DETAIL_TASK_QUEUE_H_

#include <optional>
#include <tuple>
#include <variant>

#include "window.h"

namespace camel::detail {

using TaskIdType = std::uint64_t;

using AlignmentResult = EdlibAlignResult;
using WindowConstructResult = std::vector<ReferenceWindow>;
using WindowConsensusResult = std::string;

using ResultVariant =
    std::variant<AlignmentResult, WindowConstructResult, WindowConsensusResult>;

struct TaskResult {
  TaskIdType task_id;
  ResultVariant value;
};

using AlignmentArgPack =
    std::tuple<std::span<std::unique_ptr<biosoup::NucleicAcid> const>,
               biosoup::Overlap>;
using WindowConstructArgPack =
    std::tuple<std::span<std::unique_ptr<biosoup::NucleicAcid> const>,
               std::span<biosoup::Overlap const>,
               std::span<EdlibAlignResult const>, std::uint32_t const,
               std::uint32_t const>;

using WindowConsensusArgPack =
    std::tuple<std::string_view, ReferenceWindowView, POAConfig>;

using ArgsPack = std::variant<AlignmentArgPack, WindowConstructArgPack,
                              WindowConsensusArgPack>;

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
