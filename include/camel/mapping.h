#ifndef CAMEL_MAPPING_H_
#define CAMEL_MAPPING_H_

#include <memory>
#include <utility>
#include <vector>

#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"
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

CAMEL_EXPORT struct OverlapInterval {
  std::uint32_t target_id;
  std::size_t first_index;
  std::size_t last_index;
};

CAMEL_EXPORT struct ReadOverlaps {
  std::vector<biosoup::Overlap> target_overlaps;
  std::vector<OverlapInterval> remote_intervals;
};

/**
 * @brief Find overlaps between reads
 */
CAMEL_EXPORT [[nodiscard]] auto FindOverlaps(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool, MapCfg const map_cfg,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads)
    -> std::vector<ReadOverlaps>;

}  // namespace camel

#endif /* CAMEL_MAPPING_H_ */
