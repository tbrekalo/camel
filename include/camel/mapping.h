#ifndef CAMEL_MAPPING_H_
#define CAMEL_MAPPING_H_


#include <utility>
#include <memory>
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

/**
 * @brief helper type for @ref Group and @ref FindOverlaps
 */
CAMEL_EXPORT struct ReadIdOvlpCnt {
  std::uint32_t read_id;
  std::size_t n_overlaps;
};

/**
 * @brief Intended only to be used as a return type
 *          from @ref FindOverlaps
 */
CAMEL_EXPORT struct Group {
  /**
   * @brief assumes that the read_ids are unique 
   *          and sorted in ascending order
   */
  std::vector<ReadIdOvlpCnt> read_n_ovlps; 
  std::vector<std::vector<biosoup::Overlap>> ovlp_vecs;
};

/**
 * @brief Find overlaps between reads 
 */
CAMEL_EXPORT [[nodiscard]] auto FindOverlaps(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool, MapCfg const map_cfg,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads)
    -> std::vector<Group>;

}  // namespace camel

#endif /* CAMEL_MAPPING_H_ */
