#include "camel/coverage.h"

#include "camel/mapping.h"

namespace camel {

auto CalculateCoverage(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads)
  -> std::vector<std::vector<Coverage>> {
  
}

}
