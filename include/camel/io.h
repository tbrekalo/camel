#ifndef CAMEL_IO_H_
#define CAMEL_IO_H_

#include <filesystem>
#include <memory>
#include <vector>

#include "biosoup/nucleic_acid.hpp"
#include "thread_pool/thread_pool.hpp"

/* export header */
#include "camel/export.h"
/* export header */

namespace camel {

[[nodiscard]] CAMEL_EXPORT auto LoadSequences(std::filesystem::path const& path)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>>;

[[nodiscard]] CAMEL_EXPORT auto LoadSequences(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::vector<std::filesystem::path> const& paths)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>>;

}  // namespace camel

#endif /* CAMEL_IO_H_ */
