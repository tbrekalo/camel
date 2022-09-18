#include "call_sites.h"

#include <cmath>

namespace camel::detail {

auto IsStableSite(CoverageSignals const& covg_signals,
                  std::uint32_t const covg_estimate,
                  std::uint8_t const base_code, double const alfa,
                  double const beta) -> bool {
  auto const covg = static_cast<double>(
      std::accumulate(covg_signals.val.cbegin(), covg_signals.val.cend(),
                      std::uint16_t(0), std::plus<std::uint16_t>{}) -
      covg_signals.val[CoverageSignals::kInsIdx]);

  auto const match_rate =
      static_cast<double>(covg_signals.val[base_code]) / covg;
  auto const ins_rate =
      static_cast<double>(covg_signals.val[CoverageSignals::kInsIdx]) / covg;

  return covg > covg_estimate && match_rate > alfa && ins_rate < beta;
}

auto IsSnpSite(CoverageSignals const& covg_signals,
               std::uint32_t const covg_estimate, double const dominance_ratio)
    -> bool {
  auto const covg = static_cast<double>(
      std::accumulate(covg_signals.val.cbegin(), covg_signals.val.cend() + 4U,
                      std::uint16_t(0), std::plus<std::uint16_t>{}));

  for (auto i = 0U; i < 4U; ++i) {
    auto const ratio = static_cast<double>(covg_signals.val[i]) / covg;
    if (ratio > dominance_ratio) {
      return false;
    }
  }

  return true;
}

}  // namespace camel::detail
