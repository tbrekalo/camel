#ifndef CAMEL_DETAIL_INTERVAL_H_
#define CAMEL_DETAIL_INTERVAL_H_

#include <cstdint>

struct Interval {
  std::uint32_t first;
  std::uint32_t last;
  friend constexpr auto operator==(Interval const lhs,
                                   Interval const rhs) noexcept -> bool {
    return lhs.first == rhs.first && lhs.last == rhs.last;
  }

  friend constexpr auto operator!=(Interval const lhs, Interval const rhs)
      -> bool {
    return !(lhs == rhs);
  }
};

constexpr auto IntervalLength(Interval const intv) -> std::uint32_t {
  return intv.last - intv.first;
}

constexpr auto LocalizeInterval(std::uint32_t const pos, Interval const intv)
    -> Interval {
  return {.first = intv.first - pos, .last = intv.last - pos};
}

#endif /* CAMEL_DETAIL_INTERVAL_H_ */
