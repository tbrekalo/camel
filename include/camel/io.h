#ifndef CAMEL_IO_H_
#define CAMEL_IO_H_

#include <filesystem>
#include <memory>
#include <vector>

#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"

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
      std::vector<std::filesystem::path> const& paths)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>>;

/**
 * @brief multithreaded sequence storage
 */
CAMEL_EXPORT auto StoreSequences(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& seqs,
    std::filesystem::path const& dst_file) -> void;

/**
 * @brief loads overlaps, keeps the best one for each read
 */
CAMEL_EXPORT [[nodiscard]] auto LoadOverlaps(
    std::filesystem::path const& paf_path,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads,
    std::size_t const n_overlaps,
    std::size_t const min_overlap_len)
    -> std::vector<std::vector<biosoup::Overlap>>;

}  // namespace camel

#endif /* CAMEL_IO_H_ */
