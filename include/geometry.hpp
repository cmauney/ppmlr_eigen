#pragma once

#include <Eigen/Dense>

#include "global_const.hpp"
#include "global_def.hpp"

namespace geo {
// geometry IDs
enum GEO { CART = 0, CYLIN, SPHERI, CYLIN_THETA, SPHERI_THETA, SPHERI_PHI };

struct planar_geo_t {
  constexpr planar_geo_t() {}

  template <typename V>
  constexpr inline V volume(const V& /*x*/, const V& dx) const {
    return dx;
  }

  template <typename V>
  constexpr inline V area_mid(const V& /*x0*/, const V& /*xn*/) const {
    return V::Ones();
  }

  template <typename V>
  constexpr inline V overlap_volume(const V& x0, const V& xn) const {
    return xn - x0;
  }
};
/*
struct cylin_geo_t {
  constexpr cylin_geo_t() {}
  constexpr inline double volume(const double x, const double dx) const {
    return dx * (x + half * dx);
  }
  constexpr inline double area_mid(const double x0, const double xn) const {
    return half * (x0 + xn);
  }
  constexpr inline double overlap_volume(const double x0,
                                         const double xn) const {
    auto dx = xn - x0;

    return dx * (x0 + half * dx);
  }
};
*/
}  // namespace geo
