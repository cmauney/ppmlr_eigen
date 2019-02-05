#pragma once

#include <cmath>
#include <iostream>
#include <valarray>

#include "global_const.hpp"
#include "global_def.hpp"
#include "ppmlr/ppm_constant.hpp"

namespace physics {

/*struct Hydro {
  const double gamma, gammam1;
  double last_dt, dt, hdt;

  double svel;

  Hydro(double gam)
      : gamma(gam), gammam1(gam - 1), last_dt(ppm::small), svel(zero) {}

  template <typename V>
  constexpr inline double ridt(const double svel, const V& xvel, const V& dx) {
    return std::max(svel, (xvel.abs() / dx).maxCoeff());
  }

  constexpr inline void set_dt(double timestep) {
    last_dt = dt;

    dt = timestep;
    hdt = half * dt;
  }

  template <typename V>
  constexpr inline void first_dt(const V& prs,
                                 const V& rho,
                                 const V& xvel,
                                 const V& dx) {
    auto ctdx = (gamma * prs / rho).sqrt() / dx;

    auto svel = ctdx.maxCoeff();

    dt = courant / ridt(svel, xvel, dx);
    hdt = half * dt;
    last_dt = dt;
  }

  template <typename V>
  constexpr inline void dt_courant(double cdx, const V& xvel, const V& dx) {
    auto dt_n = courant / ridt(cdx, xvel, dx);

    auto dtu = last_dt * 1.1;
    set_dt(std::min(dtu, dt_n));
  }
};
*/
}  // namespace physics
