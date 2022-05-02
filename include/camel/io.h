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

/**
 * @brief load sequences from single fasta/fastq file
 *
 * @return loaded fasta/fastq sequences sorted and indexed by sequence name
 *
 */
CAMEL_EXPORT [[nodiscard]] auto LoadSequences(std::filesystem::path const& path)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>>;

/**
 * @brief multithreaded sequence loading from multiple fasta/fastq files
 *
 * @return loaded fasta/fastq sequences sorted and indexed by sequence name
 *
 */
CAMEL_EXPORT [[nodiscard]] auto LoadSequences(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::vector<std::filesystem::path> const& paths)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>>;

/**
 * @brief serialize pile batch to the destination folder in a single file
 *
 * @return path to the serialized pile batch
 *
 */
CAMEL_EXPORT auto SerializePileBatch(std::vector<Pile>::const_iterator first,
                                     std::vector<Pile>::const_iterator last,
                                     std::filesystem::path const& dst_dir,
                                     std::string const& batch_name)
    -> std::filesystem::path;

/**
 * @brief serialize coverage piles to disk using default expected file size (4
 * GiB)
 *
 */
CAMEL_EXPORT auto SerializePiles(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::vector<Pile> const& piles, std::filesystem::path const& dst_dir)
    -> void;

/**
 * @brief serialize coverage piles to disk using specified expected file size
 *
 */
CAMEL_EXPORT auto SerializePiles(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::vector<Pile> const& piles, std::filesystem::path const& dst_dir,
    std::size_t const expected_file_sz) -> void;

/**
 * @brief deserialize piles from source folder
 */
CAMEL_EXPORT [[nodiscard]] auto DeserializePiles(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::filesystem::path const& src_dir) -> std::vector<Pile>;

// TODO: consider writing single threaded API

}  // namespace camel

#endif /* CAMEL_IO_H_ */
