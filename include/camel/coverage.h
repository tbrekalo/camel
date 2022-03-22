#ifndef CAMEL_COVERAGE_H_
#define CAMEL_COVERAGE_H_

#include <cstdint>
#include <vector>

#include "biosoup/nucleic_acid.hpp"
#include "camel/export.h"
#include "camel/mapping.h"
#include "thread_pool/thread_pool.hpp"

namespace camel {

struct Coverage {
  std::uint16_t ins;
  std::uint16_t del;
  std::uint16_t mat;
  std::uint16_t mis;
};

[[nodiscard]] CAMEL_EXPORT auto CalculateCoverage(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool, MapCfg const map_cfg,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& seqs)
    -> std::vector<std::vector<Coverage>>;

}  // namespace camel
#endif /* CAMEL_COVERAGE_H_ */
