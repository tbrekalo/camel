#ifndef CAMEL_MAPPING_H_
#define CAMEL_MAPPING_H_

#include <memory>
#include <utility>
#include <vector>

#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"
#include "camel/state.h"
#include "thread_pool/thread_pool.hpp"

/* export header */
#include "camel/export.h"
/* export header */

namespace camel {

CAMEL_EXPORT struct MapCfg {
  MapCfg() = default;
  MapCfg(std::uint8_t kmer_len, std::uint8_t win_len, double filter_p);

  std::uint8_t kmer_len = 15;
  std::uint8_t win_len = 5;
  double filter_p = 0.001;
};

// TODO: write overlap storage class

/**
 * @brief Find all overlaps between reads
 */
CAMEL_EXPORT [[nodiscard]] auto FindOverlaps(
    State& state, MapCfg const map_cfg,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads)
    -> std::vector<std::vector<biosoup::Overlap>>;

}  // namespace camel

#endif /* CAMEL_MAPPING_H_ */
