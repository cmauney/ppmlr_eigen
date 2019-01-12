#pragma once

#include <cstddef>
#include <iostream>

//#include "hydro_physics.hpp"
#include "geometry.hpp"
#include "global_def.hpp"
#include "ppm_constant.hpp"
#include "ppm_def.hpp"

namespace ppm
{

template <typename GeometorType, typename Frame>
struct PPMGrid {
  using frame = Frame;
  using vector_type = VectorType<frame::iN>;
  using geo_type = GeometorType;

  // grid values
  std::array<vector_type, 2> xe, xc, dx, dvol;

  // storage space for wiggle
  vector_type __xe, __xc, __dx, __dvol;

  vector_type amid, xnolap, xndiff;

  //  sweep_vector_t<NZ_SWEEP> r, u, v, w, p, e, f;
  //  sweep_vector_t<NZ_SWEEP> q;

  // geometry variable
  geo_type geometer;
  double radius;

  PPMGrid() : radius(1.0) {}

  template <index_type I>
  constexpr inline void build_volume() {
    dvol[I] = geometer.volume(xe[I], dx[I]);
  }

  constexpr inline void build_amid() { amid = geometer.area_mid(xe[0], xe[1]); }

  constexpr inline void build_overlap() {
    xnolap = geometer.overlap_volume(xe[0], xe[1]);
    xndiff = xe[1] - xe[0];
  }

  constexpr inline void store_coords() {
    __xe.swap(xe[0]);
    __xc.swap(xc[0]);
    __dx.swap(dx[0]);
    __dvol.swap(dvol[0]);
    /*    for (index_type i = bp::i0; i < bp::iN; ++i) {*/
    //__xe[i] = xe[0][i];
    //__xc[i] = xc[0][i];
    //__dx[i] = dx[0][i];
    //__dvol[i] = dvol[0][i];
    //}

    //  __xe[vector_type::iN] = xe[0][vector_type::iN];
  }

  inline void retieve_coords() {
    xe[0].swap(__xe);
    xc[0].swap(__xe);
    dx[0].swap(__xe);
    dvol[0].swap(__xe);
    /* for (index_type i = bp::i0; i < bp::iN; ++i) {*/
    // xe[0][i] = __xe[i];
    // xc[0][i] = __xc[i];
    // dx[0][i] = __dx[i];
    // dvol[0][i] = __dvol[i];
    /*}*/
    //    xe[0][vector_type::iN] = __xe[vector_type::iN];
  }
};

} // namespace ppm
