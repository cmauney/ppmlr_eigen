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
  using geo_type = GeometorType;

  using frame::iN;
  // grid values
//  MatrixType<iN, 2> xe, xc, dx, dvol;

  VectorType<iN> xe0, xc0, dx0, dvol0;
  VectorType<iN> xen, xcn, dxn, dvoln;
  // storage space for wiggle
  VectorType<iN> __xe, __xc, __dx, __dvol;

  VectorType<iN> amid, xnolap, xndiff;

  //  sweep_vector_t<NZ_SWEEP> r, u, v, w, p, e, f;
  //  sweep_vector_t<NZ_SWEEP> q;

  // geometry variable
  geo_type geometer;
  double radius;

  PPMGrid() : radius(1.0) {}

  template<typename GridSlice>
  constexpr inline void set_t0( const GridSlice& in_g )
  {
    xe0 = in_g.pad({NZ_GHOST, NZ_GHOST});
    xc0 = in_g.pad({NZ_GHOST, NZ_GHOST});
    dx0 = in_g.pad({NZ_GHOST, NZ_GHOST});
  }

  template <index_type I>
  constexpr inline void build_volume() {
    if constexpr( I == 0 )
      dvol0 = geometer.volume(xe0, dx0);
    else
      dvoln = geometer.volume(xen, dxn);
  }

  constexpr inline void build_amid() { amid = geometer.area_mid(xe0, xen); }

  constexpr inline void build_overlap() {
    xnolap = geometer.overlap_volume(xe0, xen);
    xndiff = xen - xe0;
  }

  constexpr inline void store_coords() {
    __xe.swap(xe0);
    __xc.swap(xc0);
    __dx.swap(dx0);
    __dvol.swap(dvol0);
    /*    for (index_type i = bp::i0; i < bp::iN; ++i) {*/
    //__xe[i] = xe[0][i];
    //__xc[i] = xc[0][i];
    //__dx[i] = dx[0][i];
    //__dvol[i] = dvol[0][i];
    //}

    //  __xe[vector_type::iN] = xe[0][vector_type::iN];
  }

  inline void retieve_coords() {
    xe0.swap(__xe);
    xc0.swap(__xe);
    dx0.swap(__xe);
    dvol0.swap(__xe);
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
