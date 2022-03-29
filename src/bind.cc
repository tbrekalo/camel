#include "camel/coverage.h"
#include "camel/io.h"
#include "nanobind/nanobind.h"
#include "nanobind/stl/string.h"
#include "nanobind/stl/vector.h"

namespace nb = nanobind;

NB_MODULE(camelpy_ext, m) {
  m.def("hello", []() -> std::string { return "hello"; });
}
