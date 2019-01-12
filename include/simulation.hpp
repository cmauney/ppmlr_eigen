#pragma once

#include "ppmlr/ppm_grid.hpp"

#include "geometry.hpp"
#include "global_const.hpp"
#include "global_def.hpp"

namespace simulation {

constexpr index_type NDIM = 1;

// X dim definitions
using XGEO = geo::planar_geo_t;
constexpr index_type NX = 128;

}  // namespace simulation
