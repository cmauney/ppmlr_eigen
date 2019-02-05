#pragma once

#include <cmath>
#include <iostream>
#include <valarray>

#include "physics/eos.hpp"
#include "global_const.hpp"
#include "global_def.hpp"
//#include "ppmlr/ppm_constant.hpp"

namespace physics {

constexpr double courant = 0.5;

struct Timestep {
  double last_dt, dt, hdt;

  Timestep() {}

  template <typename V>
  constexpr inline double ridt(const double svel, const V& xvel, const V& dx) {
    return std::max(svel, (xvel.abs() / dx).maxCoeff());
  }

  constexpr inline void set_dt(double timestep) {
    last_dt = dt;

    dt = timestep;
    hdt = half * dt;
  }

  template <class V>
  constexpr inline void first_dt(const physics::IdealGas& eos,
                                 const V& prs,
                                 const V& rho,
                                 const V& xvel,
                                 const V& dx) {

    auto svel = (eos.c(prs, rho) / dx).maxCoeff();

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

}  // namespace physics
