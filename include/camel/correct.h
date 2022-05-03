#ifndef CAMEL_CORRECT_H_
#define CAMEL_CORRECT_H_

#include <memory>
#include <vector>

#include "biosoup/nucleic_acid.hpp"
#include "camel/export.h"
#include "camel/mapping.h"
#include "thread_pool/thread_pool.hpp"

namespace camel {

struct AnnotatedRead {
  std::unique_ptr<biosoup::NucleicAcid> read;
  std::vector<std::uint32_t> snp_calls;
};

CAMEL_EXPORT [[nodiscard]] auto SnpErrorCorrect(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> src_reads)
    -> std::vector<AnnotatedRead>;

}  // namespace camel

#endif /* CAMEL_CORRECT_H_ */
