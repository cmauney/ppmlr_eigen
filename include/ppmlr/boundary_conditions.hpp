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
#include <Eigen/Dense>

namespace ppm {

template<class Grid, index_type N>
struct BoundaryCondition
{
  using grid_type = Grid;
  using vector_type = __PPMVector<N>;
  using fluid_field = __PPMFluidVars<N>;
  //  using frame = typename grid_type::frame;
  //  using vector_type = VectorType<frame::iN>;

  BoundaryCondition() = default;
  ~BoundaryCondition() = default;

  constexpr inline void reflect_var(vector_type& v, double symm = 1.0)
  {
    for (index_type i = 0; i < NZ_GHOST; ++i) {
      v[NZ_GHOST - i - 1] = symm * v[NZ_GHOST + i];
      v[N + NZ_GHOST + i] = symm * v[N + NZ_GHOST - i - 1];
    }
  }

  constexpr inline void extend_grid(grid_type& g)
  {
    for (index_type i = 0; i < NZ_GHOST; ++i) {
      index_type l = NZ_GHOST - i;
      index_type r = N + NZ_GHOST + i;

      g.xe0[l - 1] = g.xe0[l] - g.dx0[l - 1];

      g.xc0[l - 1] = g.xc0[l] - g.dx0[l - 1];

      g.xe0[r] = g.xe0[r - 1] + g.dx0[r - 1];

      g.xc0[r] = g.xc0[r - 1] + g.dx0[r - 1];
    }
  }

  constexpr inline void apply(grid_type& g, fluid_field& vars)
  {
    reflect_var(vars[RHO]);
    reflect_var(vars[PRS]);
    reflect_var(vars[XVL], -1.0);
    reflect_var(vars[YVL]);
    reflect_var(vars[ZVL]);
    reflect_var(vars[ENT]);
    reflect_var(vars[FLT]);

    reflect_var(g.dx0);
    extend_grid(g);
  }
};

} // namespace ppm
