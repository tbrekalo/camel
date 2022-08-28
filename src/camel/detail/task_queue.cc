#include "task_queue.h"

#include <atomic>
#include <functional>

#include "tbb/concurrent_queue.h"
#include "tbb/task_group.h"

namespace camel::detail {

template <class... Ts>
struct overload : Ts... {
  using Ts::operator()...;
};

template <class... Ts>
overload(Ts...) -> overload<Ts...>;

struct TaskQueue::Impl {
  using Work = std::function<void()>;

  std::atomic<std::size_t> task_counter;
  tbb::task_group task_group;
  tbb::concurrent_queue<std::function<void()>> task_queue[3];
  tbb::concurrent_queue<TaskResult> result_queue;

  ~Impl() { task_group.wait(); }
};

template <class ArgsTuple>
struct ArgPackTaskTraits;

template <>
struct ArgPackTaskTraits<AlignmentArgPack> {
  static constexpr auto function = AlignStrings;
  static constexpr auto priority = 2;
};

template <>
struct ArgPackTaskTraits<WindowArgPack> {
  static constexpr auto function = CreateWindowsFromAlignments;
  static constexpr auto priority = 1;
};

TaskQueue::TaskQueue() : pimpl_(std::make_shared<Impl>()) {}

auto TaskQueue::Push(ArgsPack args_pack) -> std::size_t {
  auto task_id = pimpl_->task_counter++;
  pimpl_->task_group.run([=]() -> void {
    std::visit(
        [=, pimpl = pimpl_](auto args_pack) -> void {
          using traits = ArgPackTaskTraits<
              std::remove_cv_t<std::remove_reference_t<decltype(args_pack)>>>;
          pimpl->task_queue[traits::priority].emplace([=]() -> void {
            pimpl->result_queue.push(
                TaskResult{.task_id = task_id,
                           .value = std::apply(traits::function, args_pack)});
          }

          );
        },
        args_pack);

    {
      auto work = Impl::Work();
      for (auto idx = 0U; true; idx = (idx + 1U) % 3U) {
        if (pimpl_->task_queue[idx].try_pop(work)) {
          break;
        }
      }

      work();
    }
  });

  return task_id;
}

auto TaskQueue::TryPop() -> std::optional<TaskResult> {
  auto task_result = TaskResult();
  if (pimpl_->result_queue.try_pop(task_result)) {
    return task_result;
  }

  return std::nullopt;
}

}  // namespace camel::detail
