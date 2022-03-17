#ifndef CAMEL_CAMEL_H_
#define CAMEL_CAMEL_H_

#include <cstdint>
#include <vector>

#include "camel/export.h"

namespace camel {

struct Coverage {
  std::uint16_t ins;
  std::uint16_t del;
  std::uint16_t mat;
  std::uint16_t mis;
};



}  // namespace camel

#endif /* CAMEL_CAMEL_H_ */
