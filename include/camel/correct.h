#ifndef CAMEL_CORRECT_H_
#define CAMEL_CORRECT_H_

#include <memory>
#include <vector>

#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"
#include "camel/poa_config.h"
#include "camel/export.h"

namespace camel {

CAMEL_EXPORT struct CorrectConfig {
  POAConfig poa_cfg;
  std::uint32_t window_length;
};

CAMEL_EXPORT [[nodiscard]] auto ErrorCorrect(
    CorrectConfig const correct_cfg,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> src_reads,
    std::vector<std::vector<biosoup::Overlap>> overlaps)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>>;

}  // namespace camel

#endif /* CAMEL_CORRECT_H_ */
