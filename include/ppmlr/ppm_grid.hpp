#pragma once

#include <cstddef>
#include <iostream>

//#include "hydro_physics.hpp"
#include "geometry.hpp"
#include "global_def.hpp"
#include "ppm_constant.hpp"
#include "ppm_def.hpp"

namespace ppm {

template<class GeometorType, index_type N>
struct PPMGrid
{
  //  using frame = Frame;
  using geo_type = GeometorType;
  using vector_type = __PPMVector<N>;
  // grid values
  vector_type xe0, xc0, dx0, dvol0;
  vector_type xen, xcn, dxn, dvoln;

  // storage space for wiggle
  vector_type __xe, __xc, __dx, __dvol;

  vector_type amid, xnolap, xndiff;

  // geometry variable
  geo_type geometer;
  double radius;

  PPMGrid()
    : radius(1.0)
  {}

  template<class V>
  constexpr inline void set_g0(const V& xe, const V& xc, const V& dx)
  {
    auto r1 = Eigen::seqN(NZ_GHOST, Eigen::fix<N>);
    xe0(r1) = xe;
    xc0(r1) = xc;
    dx0(r1) = dx;
  }

  template<index_type I>
  constexpr inline void build_volume()
  {
    if constexpr (I == 0)
      dvol0 = geometer.volume(xe0, dx0);
    else
      dvoln = geometer.volume(xen, dxn);
  }

  constexpr inline void build_amid() { amid = geometer.area_mid(xe0, xen); }

  constexpr inline void build_overlap()
  {
    xnolap = geometer.overlap_volume(xe0, xen);
    xndiff = xen - xe0;
  }

  constexpr inline void store_coords()
  {
    __xe.swap(xe0);
    __xc.swap(xc0);
    __dx.swap(dx0);
    __dvol.swap(dvol0);
  }

  constexpr inline void retieve_coords()
  {
    xe0.swap(__xe);
    xc0.swap(__xe);
    dx0.swap(__xe);
    dvol0.swap(__xe);
  }
};

} // namespace ppm
