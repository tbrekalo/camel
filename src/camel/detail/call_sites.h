#ifndef CAMEL_DETAIL_CALL_SITES_H_
#define CAMEL_DETAIL_CALL_SITES_H_

#include "coverage.h"

namespace camel::detail {

[[nodiscard]] auto IsUnstableSite(CoverageSignals const& covg,
                                  std::uint32_t const covg_estimate,
                                  double const strong_base_ratio,
                                  double const indle_ratio) -> bool;

[[nodiscard]] auto IsSnpSite(CoverageSignals const& covg,
                             std::uint32_t const covg_estimate) -> bool;

[[nodiscard]] auto IsErrorSite(CoverageSignals const& covg,
                               std::uint32_t const covg_estimate) -> bool;

template <class PredFn>
[[nodiscard]] auto CallSites(std::vector<CoverageSignals> const& covg,
                             std::uint32_t const covg_estimate, PredFn pred_fn)
    -> std::enable_if_t<
        std::is_invocable_r_v<bool, PredFn, CoverageSignals const,
                              std::uint32_t const>,
        std::vector<std::uint32_t>> {
  auto buff = std::vector<std::uint32_t>();
  buff.reserve(buff.size() / 100);
  for (auto i = 0U; i < covg.size(); ++i) {
    if (pred_fn(covg[i], covg_estimate)) {
      buff.push_back(i);
    }
  }

  return decltype(buff)(buff.begin(), buff.end());
}

[[nodiscard]] auto CallSnpCandidates(std::vector<CoverageSignals> const& covg,
                                     std::uint32_t const covg_estimate)
    -> std::vector<std::uint32_t>;

[[nodiscard]] auto CallErrorSites(std::vector<CoverageSignals> const& covg,
                                  std::uint32_t const covg_estimate)
    -> std::vector<std::uint32_t>;

}  // namespace camel::detail

#endif /* CAMEL_DETAIL_CALL_SITES_H_ */
