#pragma once

#include <tuple>
#include <utility>

#include "global_const.hpp"
#include "global_def.hpp"

namespace physics {

struct IdealGas
{
  double gamma, gammam1;

  IdealGas(double _gamma)
    : gamma(_gamma)
    , gammam1(gamma - one)
  {}
  ~IdealGas() = default;

  template<typename V>
  constexpr inline decltype(auto) e_internal(const V& p, const V& r) const
  {
    return p / (r * gammam1);
  }

  template<typename V>
  constexpr inline decltype(auto) pressure(const V& ei, const V& r) const
  {
    return gammam1 * r * ei;
  }

  template<typename V>
  constexpr inline decltype(auto) c(const V& p, const V& r) const
  {
    return (gamma * p / r).sqrt();
  }
};

} // namespace physics
