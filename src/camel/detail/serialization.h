#ifndef CAMEL_SERIALIZATION_H_
#define CAMEL_SERIALIZATION_H_

#include <algorithm>
#include <cstring>
#include <filesystem>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

#include "camel/meta.h"

namespace camel::detail {

template <class T>
struct FatPointer {
  std::size_t n;
  T* ptr;

  auto begin() noexcept -> T* { return ptr; }
  auto begin() const noexcept -> T const* { return ptr; }

  auto end() noexcept -> T* { return ptr + n; }
  auto end() const noexcept -> T const* { return ptr + n; }
};

template <class T>
class IsTriviallyArchivable {
 private:
  using U = std::remove_cv_t<std::remove_reference_t<T>>;

 public:
  static constexpr auto value =
      std::is_object_v<U> && std::is_trivially_constructible_v<U>;
};

template <class T>
constexpr auto IsTriviallyArchivableV = IsTriviallyArchivable<T>::value;

template <class T, class U, class = std::void_t<>>
struct CanStoreOn : IsTriviallyArchivable<T> {};

template <class T, class U>
struct CanStoreOn<
    T, U,
    std::void_t<decltype(Store(std::declval<U&>(), std::declval<T const&>()))>>
    : std::true_type {};

template <class T, class U>
constexpr auto CanStoreOnV = CanStoreOn<T, U>::value;

template <class T, class U, class = std::void_t<>>
struct CanLoadFrom : IsTriviallyArchivable<T> {};

template <class T, class U>
struct CanLoadFrom<
    T, U, std::void_t<decltype(Load(std::declval<U&>(), std::declval<T&>()))>>
    : std::true_type {};

template <class T, class U>
constexpr auto CanLoadFromV = CanLoadFrom<T, U>::value;

class BinaryInBuffer;

class BinaryOutBuffer {
 public:
  BinaryOutBuffer() = default;
  BinaryOutBuffer(std::size_t capacity);

  BinaryOutBuffer(BinaryOutBuffer const&) = delete;
  BinaryOutBuffer& operator=(BinaryOutBuffer const&) = delete;

  BinaryOutBuffer(BinaryOutBuffer&&) = default;
  BinaryOutBuffer& operator=(BinaryOutBuffer&&) = default;

  ~BinaryOutBuffer() = default;

  template <class... Args>
  auto operator()(Args&&... args) -> void {
    (InvokeStore(std::forward<Args>(args)), ...);
  }

  auto Bytes() const -> std::vector<std::byte>;

 private:
  friend class BinaryInBuffer;

  template <class T>
  auto InvokeStore(T&& t) -> void {
    static_assert(
        CanStoreOnV<T, BinaryOutBuffer>,
        "T must satisfy IsTriviallyArchivable or implement Store function");
    if constexpr (IsTriviallyArchivableV<T>) {
      StoreBytes(BitCast<std::byte const*>(std::addressof(t)), sizeof(t));
    } else {
      Store(*this, std::forward<T>(t));
    }
  }

  template <class T>
  auto InvokeStore(FatPointer<T> const src_fat_ptr)
      -> std::enable_if_t<IsTriviallyArchivableV<T>> {
    StoreBytes(BitCast<std::byte const*>(std::addressof(src_fat_ptr.n)),
               sizeof(src_fat_ptr.n));
    StoreBytes(BitCast<std::byte const*>(src_fat_ptr.ptr),
               sizeof(T) * src_fat_ptr.n);
  }

  auto StoreBytes(std::byte const* src, std::size_t const n_bytes) -> void;

  std::vector<std::byte> buffer_;
};

class BinaryInBuffer {
 public:
  BinaryInBuffer(std::vector<std::byte> bytes);
  BinaryInBuffer(BinaryOutBuffer& binary_out);

  BinaryInBuffer(BinaryInBuffer const&) = delete;
  BinaryInBuffer& operator=(BinaryInBuffer const&) = delete;

  BinaryInBuffer(BinaryInBuffer&&) = default;
  BinaryInBuffer& operator=(BinaryInBuffer&&) = default;

  ~BinaryInBuffer() = default;

  template <class... Args>
  auto operator()(Args&&... args) -> void {
    static_assert(std::conjunction_v<std::is_lvalue_reference<Args>...,
                                     std::negation<std::is_const<Args>>...>,
                  "All parameters must be non const lvalue references");
    (InvokeLoad(std::forward<Args>(args)), ...);
  }

 private:
  template <class T>
  auto InvokeLoad(T& t) -> void {
    static_assert(
        CanLoadFromV<T, BinaryInBuffer>,
        "T must satisfy IsTriviallyArchivable or implement Load function");
    if constexpr (IsTriviallyArchivableV<T>) {
      LoadBytes(BitCast<std::byte*>(std::addressof(t)), sizeof(T));
    } else {
      Load(*this, t);
    }
  }

  template <class T>
  auto InvokeLoad(FatPointer<T>& dst_fat_ptr)
      -> std::enable_if_t<IsTriviallyArchivableV<T>> {
    LoadBytes(BitCast<std::byte*>(std::addressof(dst_fat_ptr.n)),
              sizeof(dst_fat_ptr.n));

    dst_fat_ptr.ptr = BitCast<T*>(std::addressof(buffer_[pos_]));
    ShiftPos(sizeof(T) * dst_fat_ptr.n);
  }

  auto LoadBytes(std::byte* dst, std::size_t const n_bytes) -> void;
  auto ShiftPos(std::size_t const n_bytes) -> void;

  std::vector<std::byte> buffer_;
  std::size_t pos_;
};

template <class Buff>
auto Store(Buff& buff, std::string const& str) -> void {
  buff(FatPointer<std::add_const_t<std::string::value_type>>{
      .n = str.size(), .ptr = str.data()});
}

template <class Buff>
auto Load(Buff& buff, std::string& dst) -> void {
  auto buff_view = FatPointer<std::string::value_type>();
  buff(buff_view);

  dst.clear();
  dst.resize(buff_view.n);
  std::copy(buff_view.begin(), buff_view.end(), dst.begin());
}

template <class Buff, class T>
auto Store(Buff& buff, std::vector<T> const& vec)
    -> std::enable_if_t<IsTriviallyArchivableV<T>> {
  buff(FatPointer<std::add_const_t<T>>{.n = vec.size(), .ptr = vec.data()});
}

template <class Buff, class T>
auto Load(Buff& buff, std::vector<T>& dst)
    -> std::enable_if_t<IsTriviallyArchivableV<T>> {
  auto buff_view = FatPointer<T>{};
  buff(buff_view);

  dst.clear();
  dst.resize(buff_view.n);
  std::copy(buff_view.begin(), buff_view.end(), dst.begin());
}

template <class Buff, class T>
auto Store(Buff& buff, std::vector<T> const& vec)
    -> std::enable_if_t<!IsTriviallyArchivableV<T>> {
  buff(vec.size());
  for (auto const& it : vec) {
    buff(it);
  }
}

template <class Buff, class T>
auto Load(Buff& buff, std::vector<T>& dst)
    -> std::enable_if_t<!IsTriviallyArchivableV<T>> {
  auto n = std::size_t();
  buff(n);

  dst.clear();
  dst.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    buff(dst[i]);
  }
}

auto CompressToFile(std::vector<std::byte> const& bytes,
                    std::filesystem::path const& path) -> void;

[[nodiscard]] auto DecompressFromFile(std::filesystem::path const& path)
    -> std::vector<std::byte>;

}  // namespace camel::detail

#endif /* CAMEL_SERIALIZATION_H_ */
