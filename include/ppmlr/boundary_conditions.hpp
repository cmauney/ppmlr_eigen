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

template <class Grid>
struct BoundaryCondition {
  using grid_type = Grid;
  using frame = typename grid_type::frame;
//  using vector_type = VectorType<frame::iN>;

  BoundaryCondition() = default;
  ~BoundaryCondition() = default;

  template<class V>
  void reflect_bc(V& v, index_type c, double symm = 1.0) {
    for (index_type i = frame::i0; i < frame::j0; ++i) {
      v.col(c)[frame::j0 - i - 1] = symm * v.col(c)[i + frame::j0];
    }

    for (index_type i = frame::i0; i < frame::j0; ++i) {
      v.col(c)[frame::jM + i] = symm * v.col(c)[frame::jM - i - 1];
    }
  }

  template<class V>
  inline void apply(grid_type& g, V& vars) {
    reflect_bc(vars,RHO);
    reflect_bc(vars,PRS);
    reflect_bc(vars,XVL, -1.0);
    reflect_bc(vars,YVL);
    reflect_bc(vars,ZVL);
    reflect_bc(vars,ENT);
    reflect_bc(vars,FLT);

    for (index_type i = frame::i0; i < frame::j0; ++i) {
      g.dx0[frame::j0-i-1] = g.dx0[i+frame::j0];
      g.dx0[frame::jM+i] = g.dx0[frame::jM+i-1];

      g.xe0[frame::j0 - i - 1] =
          g.xe0[frame::j0 - i] - g.dx0[frame::j0 - i - 1];

      g.xc0[frame::j0 - i - 1] =
          g.xc0[frame::j0 - i] - g.dx0[frame::j0 - i - 1];

      g.xe0[frame::jM + i] =
          g.xe0[frame::jM + i - 1] + g.dx0[frame::jM + i - 1];
      
      g.xc0[frame::jM + i] =
          g.xc0[frame::jM + i - 1] + g.dx0[frame::jM + i - 1];


    }
  }
};

}  // namespace ppm
