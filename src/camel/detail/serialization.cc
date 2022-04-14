#include "serialization.h"

#include <fstream>
#include <stdexcept>

#include "camel/meta.h"
#include "fmt/core.h"
#include "zstd.h"

namespace camel::detail {

constexpr auto kLibBufferCapacity = 1U << 22U;

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

auto CompressToFile(std::vector<std::byte> const& bytes,
                    std::filesystem::path const& path) -> void {
  auto dst_fstrm =
      std::fstream(path, std::ios::out | std::ios::binary | std::ios::trunc);

  auto compress_buff = std::vector<std::byte>(bytes.size());
  auto const compress_status =
      ZSTD_compress(compress_buff.data(), compress_buff.size(), bytes.data(),
                    bytes.size(), 5U);

  if (ZSTD_isError(compress_status)) {
    throw std::runtime_error("[camel::CompressToFile] failed to compress " +
                             path.string());
  }

  dst_fstrm.write(BitCast<char const*>(compress_buff.data()), compress_status);
}

auto DecompressFromFile(std::filesystem::path const& path)
    -> std::vector<std::byte> {
  auto in_buff = std::vector<std::byte>(std::filesystem::file_size(path));
  auto src_fstrm = std::fstream(path, std::ios::in | std::ios::binary);

  src_fstrm.read(BitCast<char*>(in_buff.data()), in_buff.size());

  // TODO: make more robust to errors
  auto const decom_sz =
      ZSTD_getFrameContentSize(in_buff.data(), in_buff.size());
  auto dst = std::vector<std::byte>(decom_sz);

  auto decompress_status =
      ZSTD_decompress(dst.data(), dst.size(), in_buff.data(), in_buff.size());

  if (ZSTD_isError(decompress_status)) {
    throw std::runtime_error("[camel::DecompressToFile] failed to decompress " +
                             path.string());
  }

  return dst;
}

}  // namespace camel::detail
