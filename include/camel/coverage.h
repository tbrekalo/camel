#ifndef CAMEL_COVERAGE_H_
#define CAMEL_COVERAGE_H_

#include <cstdint>
#include <vector>

#include "biosoup/nucleic_acid.hpp"
#include "camel/export.h"
#include "thread_pool/thread_pool.hpp"

namespace camel {

struct Coverage {
  std::uint16_t ins;
  std::uint16_t del;
  std::uint16_t mat;
  std::uint16_t mis;
};

CAMEL_EXPORT auto CalculateCoverage(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads)
    -> std::vector<std::vector<Coverage>>;

}  // namespace camel
#endif /* CAMEL_COVERAGE_H_ */
