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

CAMEL_EXPORT struct ReadOverlapsPair {
  std::unique_ptr<biosoup::NucleicAcid> read;
  std::vector<biosoup::Overlap> overlaps;
};

CAMEL_EXPORT struct Coverage {
  using ValueType = std::uint16_t;

  ValueType a;
  ValueType c;
  ValueType g;
  ValueType t;

  // ValueType mat;
  // ValueType mis;
  ValueType del;
  ValueType ins;
};

CAMEL_EXPORT struct Pile {
  std::uint32_t id;
  std::string seq_name;
  std::vector<Coverage> covgs;
};

CAMEL_EXPORT auto CalculateCoverage(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::vector<ReadOverlapsPair> const& reads_overlaps,
    std::filesystem::path const& pile_storage_dir) -> void;

}  // namespace camel

#endif /* CAMEL_COVERAGE_H_ */
