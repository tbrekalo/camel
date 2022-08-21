#include "call_sites.h"

#include <cmath>

namespace camel::detail {

static auto constexpr kSqrt5 = 2.2360679775;

[[nodiscard]] auto IsUnstableSite(CoverageSignals const& covg,
                                  std::uint32_t const covg_estimate,
                                  double const strong_base_ratio,
                                  double const indle_ratio) -> bool {
  auto const sig_sum =
      std::accumulate(covg.signals.cbegin(), covg.signals.cend(), 0UL,
                      std::plus<std::uint32_t>());

  if (*std::max_element(covg.signals.cbegin(), covg.signals.cend()) <
      covg_estimate / 3) {
    return false;
  }

  auto const strong_base_threshold =
      static_cast<std::uint32_t>(std::round(strong_base_ratio * sig_sum));
  auto const indle_threshold =
      static_cast<std::uint32_t>(std::round(indle_ratio * sig_sum));
  auto const cutoff_threshold =
      static_cast<std::uint32_t>(std::round(covg_estimate * kSqrt5));

  auto dst = true;
  for (auto i = 0U; dst && i < 4U; ++i) {
    if (covg.signals[i] >= strong_base_threshold) {
      dst = false;
    }
  }

  for (auto i = 4U; dst && i < 6U; ++i) {
    if (covg.signals[i] > indle_threshold) {
      dst = false;
    }
  }

  return dst;
}

auto IsSnpSite(CoverageSignals const& covg, std::uint32_t const covg_estimate)
    -> bool {
  auto signals = covg.signals;
  std::sort(signals.begin(), std::next(signals.begin(), 4),
            std::greater<std::uint16_t>());

  auto const sig_sum = std::accumulate(signals.cbegin(), signals.cend(), 0U,
                                       std::plus<std::uint32_t>());

  if (sig_sum > 0.5 * covg_estimate && signals[0] < 0.8 * sig_sum &&
      signals[1] > 0.2 * sig_sum && signals[0] + signals[1] >= 0.9 * sig_sum) {
    return true;
  }

  return false;
}

auto IsErrorSite(CoverageSignals const& covg, std::uint32_t const covg_estimate)
    -> bool {
  return IsUnstableSite(covg, covg_estimate, 0.95, 0.25);
}

auto CallSnpCandidates(std::vector<CoverageSignals> const& covg,
                       std::uint32_t const covg_estimate)
    -> std::vector<std::uint32_t> {
  return CallSites(covg, covg_estimate, IsSnpSite);
}

auto CallErrorSites(std::vector<CoverageSignals> const& covg,
                    std::uint32_t const covg_estimate)
    -> std::vector<std::uint32_t> {
  return CallSites(covg, covg_estimate, IsErrorSite);
}

}  // namespace camel::detail
