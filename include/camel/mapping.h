#ifndef CAMEL_MAPPING_H_
#define CAMEL_MAPPING_H_

#include <memory>
#include <vector>

#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"
#include "thread_pool/thread_pool.hpp"

namespace camel {

struct MapCfg {
  std::uint8_t kmer_len = 15;
  std::uint8_t win_len = 5;
  double filter_p = 0.01;
};

[[nodiscard]] auto FindOverlaps(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& seqs,
    MapCfg map_cfg) -> std::vector<std::vector<biosoup::Overlap>>;

}  // namespace camel

#endif /* CAMEL_MAPPING_H_ */
