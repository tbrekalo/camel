#include "nucleic_view.h"

namespace camel::detail {

NucleicView::NucleicView(biosoup::NucleicAcid const* nucleic_acid,
                         bool is_reverse_complement)
    : nucleic_acid_(nucleic_acid),
      fetch_code_impl_(is_reverse_complement ? &FetchReverseComplementCodeImpl
                                             : &FetchCodeImpl) {}

auto NucleicView::InflatedLenght() const noexcept -> std::uint32_t {
  return nucleic_acid_->inflated_len;
}

auto NucleicView::Code(std::size_t const pos) const noexcept -> std::uint8_t {
  return fetch_code_impl_(nucleic_acid_, pos);
}

auto NucleicView::InflateData(std::uint32_t const pos,
                              std::uint32_t const len) const -> std::string {
  auto dst = std::string(len, '\0');
  for (auto i = 0U; i < len; ++i) {
    dst[i] = biosoup::kNucleotideDecoder[Code(pos + i)];
  }
  return dst;
}

auto NucleicView::FetchCodeImpl(biosoup::NucleicAcid const* nucleic_acid,
                                std::size_t pos) noexcept -> std::uint8_t {
  return ((nucleic_acid->deflated_data[pos >> 5] >> ((pos << 1) & 63)) & 3);
}

auto NucleicView::FetchReverseComplementCodeImpl(
    biosoup::NucleicAcid const* nucleic_acid, std::size_t pos) noexcept
    -> std::uint8_t {
  return FetchCodeImpl(nucleic_acid, nucleic_acid->inflated_len - 1 - pos) ^ 3;
}

}  // namespace camel::detail
