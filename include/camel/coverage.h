#ifndef CAMEL_COVERAGE_H_
#define CAMEL_COVERAGE_H_

#include <cstdint>
#include <filesystem>
#include <vector>

#include "biosoup/nucleic_acid.hpp"
#include "camel/export.h"
#include "camel/mapping.h"
#include "thread_pool/thread_pool.hpp"

namespace camel {

CAMEL_EXPORT struct Coverage {
  std::uint16_t mat;
  std::uint16_t del;
  std::uint16_t ins;
  std::uint16_t mis;
};

CAMEL_EXPORT struct Pile {
  std::uint32_t id;
  std::string seq_name;
  std::vector<Coverage> covgs;
};

[[nodiscard]] CAMEL_EXPORT auto CalculateCoverage(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool, MapCfg const map_cfg,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& seqs)
    -> std::vector<Pile>;

auto CAMEL_EXPORT
SerializePiles(std::shared_ptr<thread_pool::ThreadPool> thread_pool,
               std::vector<Pile> const& piles,
               std::filesystem::path const& dst_dir) -> void;

auto CAMEL_EXPORT SerializePiles(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::vector<Pile> const& piles, std::filesystem::path const& dst_dir,
    std::size_t const expected_file_sz) -> void;

auto CAMEL_EXPORT
DeserializePiles(std::shared_ptr<thread_pool::ThreadPool> thread_pool,
                 std::filesystem::path const& src_dir) -> std::vector<Pile>;

}  // namespace camel

#endif /* CAMEL_COVERAGE_H_ */
