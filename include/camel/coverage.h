#ifndef CAMEL_COVERAGE_H_
#define CAMEL_COVERAGE_H_

#include <cstdint>
#include <filesystem>
#include <vector>

#include "biosoup/nucleic_acid.hpp"
#include "camel/mapping.h"
#include "thread_pool/thread_pool.hpp"

/* export header */
#include "camel/export.h"
/* export header */

namespace camel {

CAMEL_EXPORT struct Coverage {
  std::uint16_t mat;
  std::uint16_t del;
  std::uint16_t ins;
  std::uint16_t mis;
};

CAMEL_EXPORT struct Pile {
  std::uint32_t id;
  std::string seq_name;
  std::vector<Coverage> covgs;
};

CAMEL_EXPORT [[nodiscard]] auto CalculateCoverage(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool, MapCfg const map_cfg,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& seqs)
    -> std::vector<Pile>;

}  // namespace camel

#endif /* CAMEL_COVERAGE_H_ */
