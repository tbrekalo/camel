#ifndef CAMEL_IO_H_
#define CAMEL_IO_H_

#include <filesystem>
#include <memory>
#include <vector>

#include "biosoup/nucleic_acid.hpp"
#include "camel/coverage.h"
#include "thread_pool/thread_pool.hpp"

/* export header */
#include "camel/export.h"
/* export header */

namespace camel {

CAMEL_EXPORT [[nodiscard]] auto LoadSequences(std::filesystem::path const& path)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>>;

CAMEL_EXPORT [[nodiscard]] auto LoadSequences(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::vector<std::filesystem::path> const& paths)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>>;

CAMEL_EXPORT auto SerializePileBatch(std::vector<Pile>::const_iterator first,
                                     std::vector<Pile>::const_iterator last,
                                     std::filesystem::path const& dst_dir,
                                     std::string const& batch_name) -> std::filesystem::path;

CAMEL_EXPORT auto SerializePiles(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::vector<Pile> const& piles, std::filesystem::path const& dst_dir)
    -> void;

CAMEL_EXPORT auto SerializePiles(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::vector<Pile> const& piles, std::filesystem::path const& dst_dir,
    std::size_t const expected_file_sz) -> void;

CAMEL_EXPORT [[nodiscard]] auto DeserializePiles(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::filesystem::path const& src_dir) -> std::vector<Pile>;

}  // namespace camel

#endif /* CAMEL_IO_H_ */
