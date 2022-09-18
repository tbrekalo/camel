#ifndef CAMEL_DETAIL_CALL_SITES_H_
#define CAMEL_DETAIL_CALL_SITES_H_

#include "coverage.h"

namespace camel::detail {

auto IsStableSite(CoverageSignals const& covg_signals,
                  std::uint32_t const covg_estimate,
                  std::uint8_t const base_code, double const min_match_rate,
                  double const max_insertion_rate) -> bool;

auto IsSnpSite(CoverageSignals const& covg_signals,
               std::uint32_t const covg_estimate, double const dominance_ratio)
    -> bool;

}  // namespace camel::detail

#endif /* CAMEL_DETAIL_CALL_SITES_H_ */
