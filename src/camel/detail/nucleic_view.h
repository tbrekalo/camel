#ifndef CAMEL_DETAIL_NUCLIEC_VIEW_H_
#define CAMEL_DETAIL_NUCLIEC_VIEW_H_

#include "biosoup/nucleic_acid.hpp"

namespace camel::detail {
class NucleicView {
 public:
  NucleicView(biosoup::NucleicAcid const* nucleic_acid,
              bool is_reverse_complement);

  auto InflatedLenght() const noexcept -> std::uint32_t;

  auto Code(std::size_t const pos) const noexcept -> std::uint8_t;

  auto InflateData(std::uint32_t const pos, std::uint32_t const len) const
      -> std::string;

 private:
  using FetchCodeImplPtr = std::uint8_t (*)(biosoup::NucleicAcid const*,
                                            std::size_t);

  static auto FetchCodeImpl(biosoup::NucleicAcid const* nucleic_acid,
                            std::size_t pos) noexcept -> std::uint8_t;

  static auto FetchReverseComplementCodeImpl(
      biosoup::NucleicAcid const* nucleic_acid, std::size_t pos) noexcept
      -> std::uint8_t;

  biosoup::NucleicAcid const* nucleic_acid_;
  FetchCodeImplPtr fetch_code_impl_;
};

}  // namespace camel::detail

#endif /* CAMEL_DETAIL_NUCLIEC_VIEW_H_ */
