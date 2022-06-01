#ifndef CAMEL_CORRECT_H_
#define CAMEL_CORRECT_H_

#include <memory>
#include <vector>

#include "biosoup/nucleic_acid.hpp"
#include "camel/export.h"
#include "camel/mapping.h"
#include "camel/state.h"
#include "thread_pool/thread_pool.hpp"

namespace camel {

CAMEL_EXPORT struct POAConfig {
  std::int8_t match = 3;
  std::int8_t mismatch = -5;
  std::int8_t gap = -3;
};

CAMEL_EXPORT struct PolishConfig {
  POAConfig poa_cfg;
};

CAMEL_EXPORT [[nodiscard]] auto SnpErrorCorrect(
    State& state, MapCfg const map_cfg, PolishConfig const polish_cfg,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> src_reads)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>>;

}  // namespace camel

#endif /* CAMEL_CORRECT_H_ */
