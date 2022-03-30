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
  using ValueType = std::uint16_t;

  Coverage() = default;
  Coverage(ValueType mat, ValueType del, ValueType ins, ValueType mis)
      : mat(mat), del(del), ins(ins), mis(mis) {}

  ValueType mat;
  ValueType del;
  ValueType ins;
  ValueType mis;
};

CAMEL_EXPORT struct Pile {
  Pile() = default;
  Pile(std::uint32_t id, std::string seq_name, std::vector<Coverage> covgs)
      : id(id), seq_name(std::move(seq_name)), covgs(std::move(covgs)) {}

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
