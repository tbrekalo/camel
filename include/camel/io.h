#ifndef CAMEL_IO_H_
#define CAMEL_IO_H_

#include <filesystem>
#include <memory>
#include <vector>

#include "biosoup/nucleic_acid.hpp"
#include "camel/coverage.h"
#include "camel/state.h"
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
    State& state, std::vector<std::filesystem::path> const& paths)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>>;

/**
 * @brief multithreaded sequence storage
 */
CAMEL_EXPORT auto StoreSequences(
    State& state,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& seqs,
    std::filesystem::path const& dst_folder) -> void;
/**
 * @brief multithreaded sequence storage
 */
CAMEL_EXPORT auto StoreSequences(
  State& state,
  std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& seqs,
  std::filesystem::path const& dst_folder,
  std::uint64_t dst_file_cap) -> void;

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
CAMEL_EXPORT auto SerializePiles(State& state, std::vector<Pile> const& piles,
                                 std::filesystem::path const& dst_dir) -> void;

/**
 * @brief serialize coverage piles to disk using specified expected file size
 *
 */
CAMEL_EXPORT auto SerializePiles(State& state, std::vector<Pile> const& piles,
                                 std::filesystem::path const& dst_dir,
                                 std::size_t const expected_file_sz) -> void;

/**
 * @brief deserialize piles from source folder
 */
CAMEL_EXPORT [[nodiscard]] auto DeserializePiles(
    State& state, std::filesystem::path const& src_dir) -> std::vector<Pile>;

// TODO: consider writing single threaded API

}  // namespace camel

#endif /* CAMEL_IO_H_ */
