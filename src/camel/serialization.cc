#include "camel/serialization.h"

#include "zlib.h"

namespace camel {

namespace detail {

constexpr auto kLibBufferCapacity = 1U << 18U;

struct GzFileDeleter {
  auto operator()(gzFile file_pr) const noexcept -> void { gzclose(file_pr); }
};

using GzFileHandle = std::unique_ptr<gzFile_s, GzFileDeleter>;

}  // namespace detail

BinaryOutBuffer::BinaryOutBuffer(std::size_t capacity) {
  buffer_.reserve(capacity);
}

auto BinaryOutBuffer::StoreBytes(std::byte const* src,
                                 std::size_t const n_bytes) -> void {
  buffer_.insert(buffer_.end(), src, src + n_bytes);
}

auto BinaryOutBuffer::Bytes() const -> std::vector<std::byte> {
  return buffer_;
}

BinaryInBuffer::BinaryInBuffer(std::vector<std::byte> bytes)
    : buffer_(std::move(bytes)), pos_(0) {}

BinaryInBuffer::BinaryInBuffer(BinaryOutBuffer& binary_out)
    : BinaryInBuffer(std::move(binary_out.buffer_)) {}

auto BinaryInBuffer::LoadBytes(std::byte* dst, std::size_t const n_bytes)
    -> void {
  std::memcpy(dst, std::addressof(buffer_[pos_]), n_bytes);
  ShiftPos(n_bytes);
}

auto BinaryInBuffer::ShiftPos(std::size_t const n_bytes) -> void {
  pos_ += n_bytes;
}

auto GzStoreBytes(std::vector<std::byte> const& bytes,
                  std::filesystem::path const& path) -> void {
  if (std::filesystem::exists(path.parent_path())) {
  }
  auto dst_file = detail::GzFileHandle(gzopen(path.c_str(), "wb"));
  gzwrite(dst_file.get(), bytes.data(),
          static_cast<std::uint32_t>(bytes.size()));
}

auto GzLoadBytes(std::filesystem::path const& path) -> std::vector<std::byte> {
  auto dst = std::vector<std::byte>();
  auto src_file = detail::GzFileHandle(gzopen(path.c_str(), "rb"));
  gzbuffer(src_file.get(), detail::kLibBufferCapacity);

  auto buffer = std::vector<std::byte>(detail::kLibBufferCapacity);
  while (true) {
    auto const n_read =
        gzread(src_file.get(), buffer.data(), detail::kLibBufferCapacity);

    // TODO: could be better
    std::copy(buffer.begin(), std::next(buffer.begin(), n_read),
              std::back_inserter(dst));

    if (static_cast<std::uint32_t>(n_read) < buffer.size()) {
      break;
    }
  }

  return dst;
}

}  // namespace camel
