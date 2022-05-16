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

CAMEL_EXPORT struct AnnotatedRead {
  std::unique_ptr<biosoup::NucleicAcid> read;
  std::vector<std::uint32_t> snp_calls;
};

CAMEL_EXPORT [[nodiscard]] auto SnpErrorCorrect(
    State& state, std::vector<std::unique_ptr<biosoup::NucleicAcid>> src_reads)
    -> std::vector<AnnotatedRead>;

}  // namespace camel

#endif /* CAMEL_CORRECT_H_ */