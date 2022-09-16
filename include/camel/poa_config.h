#ifndef CAMEL_POA_CONFIG_H_
#define CAMEL_POA_CONFIG_H_

#include <cstdint>

#include "camel/export.h"

namespace camel {

CAMEL_EXPORT struct POAConfig {
  std::int8_t match = 3;
  std::int8_t mismatch = -5;
  std::int8_t gap = -4;
};

}

#endif /* CAMEL_POA_CONFIG_H_ */
