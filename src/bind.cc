#include "camel/coverage.h"
#include "camel/io.h"
#include "camel/mapping.h"
#include "camel/state.h"
#include "nanobind/nanobind.h"
#include "nanobind/stl/pair.h"
#include "nanobind/stl/shared_ptr.h"
#include "nanobind/stl/string.h"
#include "nanobind/stl/tuple.h"
#include "nanobind/stl/unique_ptr.h"
#include "nanobind/stl/vector.h"
namespace nb = nanobind;

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0U};

NB_MODULE(camelpy_ext, m) {
  using namespace nb::literals;

  nb::class_<camel::State>(m, "State").def(nb::init<>());

  m.def(
      "init_state",
      [](std::uint32_t const n_threads,
         std::string const& log_path) -> camel::State {
        return camel::State{
            .thread_pool = std::make_shared<thread_pool::ThreadPool>(n_threads),
            .log_path = log_path};
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

  m.def(
      "load_sequences",
      [](camel::State& state, std::vector<std::string> const& paths_strs)
          -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
        auto paths = std::vector<std::filesystem::path>();

        for (auto const& it : paths_strs) {
          auto it_path = std::filesystem::path(std::move(it));
          if (std::filesystem::is_regular_file(it_path)) {
            paths.emplace_back(std::move(it_path));
          } else if (std::filesystem::is_directory(it_path)) {
            for (auto dir_entry : std::filesystem::recursive_directory_iterator(
                     std::move(it_path))) {
              if (std::filesystem::is_regular_file(dir_entry)) {
                paths.emplace_back(std::move(dir_entry));
              }
            }
          }
        }

        return camel::LoadSequences(state, paths);
      });

  nb::class_<camel::MapCfg>(m, "MapCfg")
      .def(nb::init<std::uint8_t, std::uint8_t, double>(),
           nb::arg("kmer_len") = 15, nb::arg("win_len") = 5,
           nb::arg("filter_p") = 0.01)
      .def_readwrite("kmer_len", &camel::MapCfg::kmer_len)
      .def_readwrite("win_len", &camel::MapCfg::win_len)
      .def_readwrite("filter_p", &camel::MapCfg::filter_p);

  // m.def("calculate_coverage", &camel::CalculateCoverage);

  nb::class_<camel::Coverage>(m, "Coverage")
      .def(nb::init<>())
      .def_readwrite("a", &camel::Coverage::a)
      .def_readwrite("c", &camel::Coverage::c)
      .def_readwrite("g", &camel::Coverage::g)
      .def_readwrite("t", &camel::Coverage::t)
      // .def_readwrite("match", &camel::Coverage::mat)
      // .def_readwrite("mismatch", &camel::Coverage::mis)
      .def_readwrite("deletion", &camel::Coverage::del)
      .def_readwrite("insertion", &camel::Coverage::ins);

  nb::class_<camel::Pile>(m, "Pile")
      .def(nb::init<>())
      .def_readonly("id", &camel::Pile::id)
      .def_readonly("seq_name", &camel::Pile::seq_name)
      .def_readonly("coverages", &camel::Pile::covgs);

  m.def("serialize_piles",
        [](camel::State& state, std::vector<camel::Pile> const& piles,
           std::string const& path) -> void {
          camel::SerializePiles(state, piles, std::filesystem::path(path));
        });

  m.def("deserialize_piles",
        [](camel::State& state,
           std::string const& path) -> std::vector<camel::Pile> {
          return camel::DeserializePiles(state, std::filesystem::path(path));
        });

  // utility
  m.def("sample_ins_distr",
        [](std::vector<camel::Pile> const& piles)
            -> std::vector<std::pair<std::uint16_t, std::uint64_t>> {
          auto cnt_buff = std::vector<std::uint64_t>();

          auto mx_ins = 0U;
          for (auto const& pile : piles) {
            for (auto const& covg : pile.covgs) {
              if (covg.ins > mx_ins) {
                mx_ins = covg.ins;
              }
            }
          }

          auto dst_sz = 0U;
          cnt_buff.resize(mx_ins + 1U);
          for (auto const& pile : piles) {
            for (auto const& covg : pile.covgs) {
              if (covg.ins > 0) {
                dst_sz += (++cnt_buff[covg.ins] == 1UL);
              }
            }
          }

          auto dst =
              std::vector<std::pair<std::uint16_t, std::uint64_t>>(dst_sz);
          for (auto i = 0U, j = 0U; i < cnt_buff.size(); ++i) {
            if (cnt_buff[i] > 0) {
              dst[j] = std::pair(static_cast<std::uint16_t>(i), cnt_buff[i]);
              ++j;
            }
          }

          return dst;
        });
}
