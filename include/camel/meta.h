#ifndef CAMEL_META_H_
#define CAMEL_META_H_

#include <cstring>
#include <memory>
#include <type_traits>

namespace camel {

template <class To, class From>
constexpr auto BitCast(From const& src)
    -> std::enable_if_t<sizeof(To) == sizeof(From) &&
                            std::is_trivially_copyable_v<From> &&
                            std::is_trivially_copyable_v<To>,
                        To> {
  static_assert(std::is_trivially_default_constructible_v<To>,
                "To type must be trivially default constructible");

  auto dst = To{};
  std::memcpy(std::addressof(dst), std::addressof(src), sizeof(src));

  return dst;
}

}  // namespace camel

#endif /* CAMEL_META_H_ */
