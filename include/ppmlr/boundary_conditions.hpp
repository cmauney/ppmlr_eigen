#pragma once

#include <iostream>
#include <string>
#include <tuple>
#include <utility>

//#include "gsl/gsl"

//#include "hydro_physics.hpp"
#include "ppmlr/ppm_constant.hpp"
#include "ppmlr/ppm_def.hpp"
#include "ppmlr/ppm_grid.hpp"

namespace ppm {

template <typename Grid>
struct BoundaryCondition {
  using grid_type = Grid;
  using frame = typename grid_type::frame;
  using vector_type = typename grid_type::vector_type;

  BoundaryCondition() = default;
  ~BoundaryCondition() = default;

  void reflect_bc(vector_type& v, double symm = 1.0) {
    for (index_type i = frame::i0; i < frame::j0; ++i) {
      v[frame::j0 - i - 1] = symm * v[i + frame::j0];
    }

    for (index_type i = frame::i0; i < frame::j0; ++i) {
      v[frame::jM + i] = symm * v[frame::jM - i - 1];
    }
  }

  inline void apply(grid_type& g, FluidVars<vector_type>& vars) {
    reflect_bc(g.dx[0]);
    reflect_bc(vars[RHO]);
    reflect_bc(vars[PRS]);
    reflect_bc(vars[XVL], -1.0);
    reflect_bc(vars[YVL]);
    reflect_bc(vars[ZVL]);
    reflect_bc(vars[ENT]);
    reflect_bc(vars[FLT]);

    for (index_type i = frame::i0; i < frame::j0; ++i) {
      g.xe[0][frame::j0 - i - 1] =
          g.xe[0][frame::j0 - i] - g.dx[0][frame::j0 - i - 1];

      g.xe[0][frame::jM + i] =
          g.xe[0][frame::jM + i - 1] + g.dx[0][frame::jM + i - 1];
    }
  }
};

}  // namespace ppm
