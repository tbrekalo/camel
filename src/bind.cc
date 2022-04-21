#include "camel/coverage.h"
#include "camel/io.h"
#include "camel/mapping.h"
#include "nanobind/nanobind.h"
#include "nanobind/stl/shared_ptr.h"
#include "nanobind/stl/string.h"
#include "nanobind/stl/tuple.h"
#include "nanobind/stl/unique_ptr.h"
#include "nanobind/stl/vector.h"

namespace nb = nanobind;

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0U};

NB_MODULE(camelpy_ext, m) {
  using namespace nb::literals;

  nb::class_<thread_pool::ThreadPool>(m, "ThreadPool")
      .def(nb::init<std::size_t>());

  m.def("create_thread_pool",
        [](std::size_t const n) -> std::shared_ptr<thread_pool::ThreadPool> {
          return std::make_shared<thread_pool::ThreadPool>(n);
        });

  nb::class_<biosoup::NucleicAcid>(m, "NucleicAcid")
      .def(nb::init<const std::string&, const std::string&>())
      .def(nb::init<const std::string&, const std::string&,
                    const std::string&>())
      .def(
          "__deepcopy__",
          [](const biosoup::NucleicAcid& seq,
             nb::dict) -> biosoup::NucleicAcid { return seq; },
          "memo"_a)
      .def("code", &biosoup::NucleicAcid::Code)
      .def("score", &biosoup::NucleicAcid::Score)
      .def("inflate_data", &biosoup::NucleicAcid::InflateData,
           nb::arg("i") = 0U,
           nb::arg("len") = std::numeric_limits<std::uint32_t>::max())
      .def("inflate_quality", &biosoup::NucleicAcid::InflateQuality,
           nb::arg("i") = 0U,
           nb::arg("len") = std::numeric_limits<std::uint32_t>::max())
      .def("reverse_and_complement",
           &biosoup::NucleicAcid::ReverseAndComplement)
      .def_readwrite("id", &biosoup::NucleicAcid::id)
      .def_readwrite("name", &biosoup::NucleicAcid::name)
      .def_readonly("__len__", &biosoup::NucleicAcid::inflated_len);

  m.def("set_nucleic_acid_obj_cnt", [](const std::uint32_t val) -> void {
    biosoup::NucleicAcid::num_objects = val;
  });

  m.def("load_sequences",
        [](std::shared_ptr<thread_pool::ThreadPool> thread_pool,
           std::vector<std::string> const& paths_strs)
            -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
          auto paths = std::vector<std::filesystem::path>();

          paths.reserve(paths.size());
          std::transform(paths_strs.cbegin(), paths_strs.cend(),
                         std::back_inserter(paths),
                         [](std::string const& str) -> std::filesystem::path {
                           return str;
                         });

          return camel::LoadSequences(std::move(thread_pool), paths);
        });

  nb::class_<camel::MapCfg>(m, "MapCfg")
      .def(nb::init<std::uint8_t, std::uint8_t, double>(),
           nb::arg("kmer_len") = 15, nb::arg("win_len") = 5,
           nb::arg("filter_p") = 0.01)
      .def_readwrite("kmer_len", &camel::MapCfg::kmer_len)
      .def_readwrite("win_len", &camel::MapCfg::win_len)
      .def_readwrite("filter_p", &camel::MapCfg::filter_p);

  m.def("calculate_coverage", &camel::CalculateCoverage);

  nb::class_<camel::Coverage>(m, "Coverage")
      .def(nb::init<>())
      .def_readwrite("a", &camel::Coverage::a)
      .def_readwrite("c", &camel::Coverage::c)
      .def_readwrite("g", &camel::Coverage::g)
      .def_readwrite("t", &camel::Coverage::t)
      .def_readwrite("deletion", &camel::Coverage::del)
      .def_readwrite("insertion", &camel::Coverage::ins);

  nb::class_<camel::Pile>(m, "Pile")
      .def(nb::init<>())
      .def_readonly("id", &camel::Pile::id)
      .def_readonly("seq_name", &camel::Pile::seq_name)
      .def_readonly("coverages", &camel::Pile::covgs);

  m.def("serialize_piles",
        [](std::shared_ptr<thread_pool::ThreadPool> thread_pool,
           std::vector<camel::Pile> const& piles,
           std::string const& path) -> void {
          camel::SerializePiles(std::move(thread_pool), piles,
                                std::filesystem::path(path));
        });

  m.def("deserialize_piles",
        [](std::shared_ptr<thread_pool::ThreadPool> thread_pool,
           std::string const& path) -> std::vector<camel::Pile> {
          return camel::DeserializePiles(std::move(thread_pool),
                                         std::filesystem::path(path));
        });
}
