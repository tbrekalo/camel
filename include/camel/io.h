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

CAMEL_EXPORT auto LoadSequences(std::filesystem::path const& path)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>>;

CAMEL_EXPORT auto LoadSequences(std::vector<std::filesystem::path> const& paths)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>>;

CAMEL_EXPORT auto StoreSequences(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& seqs,
    std::filesystem::path const& dst_file) -> void;

CAMEL_EXPORT auto SerializeSequences(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& seqs,
    std::ostream& ostream) -> void;

CAMEL_EXPORT auto LoadOverlaps(
    std::filesystem::path const& paf_path,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads,
    double const error_threshold, std::size_t const n_overlaps)
    -> std::vector<std::vector<biosoup::Overlap>>;

}  // namespace camel

#endif /* CAMEL_IO_H_ */
