#ifndef CAMEL_DETAIL_OVERLOAD_H_
#define CAMEL_DETAIL_OVERLOAD_H_

namespace camel::detail {

template <class... Ts>
struct overload : Ts... {
  using Ts::operator()...;
};

template <class... Ts>
overload(Ts...) -> overload<Ts...>;

}

#endif /* CAMEL_DETAIL_OVERLOAD_H_ */
