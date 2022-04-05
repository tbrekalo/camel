#include "camel/coverage.h"
#include "camel/io.h"
#include "camel/mapping.h"
#include "nanobind/nanobind.h"
#include "nanobind/stl/shared_ptr.h"
#include "nanobind/stl/string.h"
#include "nanobind/stl/tuple.h"
#include "nanobind/stl/vector.h"

namespace nb = nanobind;

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0U};

NB_MODULE(camelpy_ext, m) {
  nb::class_<camel::Coverage>(m, "Coverage")
      .def(nb::init<camel::Coverage::ValueType, camel::Coverage::ValueType,
                    camel::Coverage::ValueType,
                    camel::Coverage::ValueType>())
      .def_readwrite("match", &camel::Coverage::mat)
      .def_readwrite("deletion", &camel::Coverage::del)
      .def_readwrite("insertion", &camel::Coverage::ins)
      .def_readwrite("mismatch", &camel::Coverage::mis);

  nb::class_<camel::Pile>(m, "Pile")
      .def(nb::init<std::uint32_t, std::string const&,
                    std::vector<camel::Coverage>>())
      .def_readonly("id", &camel::Pile::id)
      .def_readonly("seq_name", &camel::Pile::seq_name)
      .def_readonly("covgs", &camel::Pile::covgs);

  nb::class_<thread_pool::ThreadPool>(m, "ThreadPool")
      .def(nb::init<std::size_t>());

  m.def("create_thread_pool",
        [](std::size_t const n) -> std::shared_ptr<thread_pool::ThreadPool> {
          return std::make_shared<thread_pool::ThreadPool>(n);
        });

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
