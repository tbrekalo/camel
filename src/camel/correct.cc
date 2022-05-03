#include "camel/correct.h"

#include <array>

#include "biosoup/timer.hpp"
#include "fmt/core.h"

namespace camel {

namespace detail {}  // namespace detail

auto SnpErrorCorrect(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> src_reads)
    -> std::vector<AnnotatedRead> {
  auto dst = std::vector<AnnotatedRead>();

  return dst;
}

}  // namespace camel
