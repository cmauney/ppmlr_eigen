#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>

#include <Eigen/Dense>

#include "global_def.hpp"
#include "physics/hydro.hpp"
#include "ppmlr/boundary_conditions.hpp"
#include "ppmlr/parabola.hpp"
#include "ppmlr/ppm_constant.hpp"
#include "ppmlr/ppm_def.hpp"
#include "ppmlr/ppm_grid.hpp"
#include "ppmlr/riemann_solve.hpp"
#include "ppmlr/state.hpp"

namespace ppm {

template <typename GeoType, typename TimestepType, typename EOSType, index_type N>
struct PPMLR {
  using frame = VectorFrame<N>;

  using grid_type = PPMGrid<GeoType, frame>;
  using bc_type = BoundaryCondition<grid_type>;
  using parabola_basis_type = ParabolaBasis<frame>;
  using parabola_type = Parabola<parabola_basis_type>;
  using state_domain_type = StateDomain<frame>;
  using state_type = State<parabola_type, state_domain_type>;

  template<index_type NP>
  using parabola_set = ParabolaSet<parabola_type, NP>;
  using state_set = StateSet<state_type>;
  using riemann = RiemannSolve<state_type>;

  using vector_type = typename grid_type::vector_type;

//  const physics::Hydro& hydro;  // reference to 'globals'
  const  TimestepType& ts;
  const EOSType& eos;

  // coordinate vars
  grid_type grid;
  // fluid vars
//  FluidVars<vector_type> fluid_v;
  PPMVarSet<N> v;
  PPMSub<N, NQV> Q;
  PPMSub<N, NUV> U;

  // boundary do-er
  bc_type bc;

  // flatten stuff
  //  flat_type flat;
  // parabolas
  parabola_set<3> parabolas_prim;
  parabola_set<6> parabolas_cons;
  // states
  state_set states;
  // riemann-er
  riemann riemann_solve;

  vector_type flat, z;
  // kinetic E store, cspeed
  vector_type ekin, cdx;
  // evolve stores
  vector_type dm, dm0, dtdm, pmid, umid, upmid, uold;
  // remap stores
  std::array<vector_type, N_FLUIDV> fluxes;
  // max svel in sweep
  double max_cdx;

  PPMLR(const TimestepType& _ts, const EOSType& _eos)
      : ts(_ts), eos(_eos), max_cdx(zero), 
      Q{v.fluid_sca.all({RHO}), v.fluid_sca.all({PRS}), v.fluid_vec.all({0, VEL})} {}

  ~PPMLR() = default;

/*  inline void put_vars(const VectorType<N>& xe,
                       const VectorType<N>& xc,
                       const VectorType<N>& dx,
                       const VectorType<N>& rho,
                       const VectorType<N>& prs,
                       const VectorType<N>& xvel,
                       const VectorType<N>& yvel,
                       const VectorType<N>& zvel,
                       const VectorType<N>& f) {*/
    template<typename GridSlice, typename ScalarSlice, typename VectorSlice>
    inline void put_vars( const GridSlice& in_g, const ScalarSlice& in_s, const VectorSlice in_v  )
    {
    auto i_lhs = Eigen::seqN(frame::j0, Eigen::fix<N>);

    set_t0(in_g);

    v.fluid_sca = in_s.pad({{NZ_GHOST, NZ_GHOST}, {0,0}});
    v.fluid_vec = in_v.pad({{NZ_GHOST, NZ_GHOST}, {0,0}, {0,0}});
    
/*    grid.xe[0](i_lhs) = in_g.xe;
    grid.xc[0](i_lhs) = in_g.xc;
    grid.dx[0](i_lhs) = in_g.dx;
*/
/*    fluid_v[RHO](i_lhs) = rho;
    fluid_v[PRS](i_lhs) = prs;
    fluid_v[XVL](i_lhs) = xvel;
    fluid_v[YVL](i_lhs) = yvel;
    fluid_v[ZVL](i_lhs) = zvel;
    fluid_v[FLT](i_lhs) = f;
*/

    build_ekin();

    fluid_v[ENT](i_lhs) =
        fluid_v[PRS](i_lhs) / (fluid_v[RHO](i_lhs) * hydro.gammam1) +
        ekin(i_lhs);
    /*for (index_type i = frame::j0; i < frame::jM; ++i)
      fluid_v[ENT][i] =
          fluid_v[PRS][i] / (fluid_v[RHO][i] * hydro.gammam1) + ekin[i];*/
  }

  inline void get_vars(VectorType<N>& rho,
                       VectorType<N>& prs,
                       VectorType<N>& xvel,
                       VectorType<N>& yvel,
                       VectorType<N>& zvel,
                       VectorType<N>& f) {
    auto i_rhs = Eigen::seqN(frame::j0, Eigen::fix<N>);
    rho = fluid_v[RHO](i_rhs);
    prs = fluid_v[PRS](i_rhs);
    xvel = fluid_v[XVL](i_rhs);
    yvel = fluid_v[YVL](i_rhs);
    zvel = fluid_v[ZVL](i_rhs);
    f = fluid_v[FLT](i_rhs);
    /*
for (index_type i = frame::j0; i < frame::jM; ++i) {
  index_type n = i - frame::j0;
  rho[n] = fluid_v[RHO][i];
  prs[n] = fluid_v[PRS][i];
  xvel[n] = fluid_v[XVL][i];
  yvel[n] = fluid_v[YVL][i];
  zvel[n] = fluid_v[ZVL][i];
  f[n] = fluid_v[FLT][i];
}*/
  }

  constexpr inline void build_ekin() {
/*    ekin = half * (fluid_v[XVL].square() + fluid_v[YVL].square() +
                   fluid_v[ZVL].square());*/
    ekin = half * (v.vector_field
  }

  constexpr inline void build_cdx() {
    auto i_ = Eigen::seqN(frame::i0 + 2, Eigen::fix<frame::iN - 4>);

    cdx(i_) = (hydro.gamma * fluid_v[PRS](i_) / fluid_v[RHO](i_)).sqrt() /
              (grid.dx[0](i_) * grid.radius);

    max_cdx = cdx(i_).maxCoeff();
  }

  inline void flatten() {
    double delp1, delp2, zsten;
    double shock, old_flat;

    for (index_type i = frame::j0 - 4; i < frame::jM + 4; ++i) {
      delp1 = fluid_v[PRS][i + 1] - fluid_v[PRS][i - 1];
      delp2 = fluid_v[PRS][i + 2] - fluid_v[PRS][i - 2];
      if (std::abs(delp2) < small)
        delp2 = small;
      shock =
          std::abs(delp1) / std::min(fluid_v[PRS][i + 1], fluid_v[PRS][i - 1]) -
          epsilon;
      shock = std::max(0.0, shock);
      if (shock > 0.0)
        shock = 1.0;
      if (fluid_v[XVL][i - 1] < fluid_v[XVL][i + 1])
        shock = 0.0;

      z[i] = shock * std::max(0.0, (delp1 / delp2 - omega1) * omega2);
    }

    z[1] = z[2];
    z[frame::iN - 2] = z[frame::iN - 3];

    flat = zero;
    //  for (const auto& [i, im2, im1, ip1, ip2] : r1) {
    for (index_type i = frame::j0 - 4; i < frame::jM + 4; ++i) {
      zsten = std::max(z[i - 1], std::max(z[i], z[i + 1]));
      flat[i] = std::max(0.0, std::min(0.5, zsten));
    }

    for (index_type i = frame::j0 - 3; i < frame::jM + 3; ++i) {
      old_flat = fluid_v[FLT][i] - static_cast<int>(fluid_v[FLT][i]);
      if (flat[i] > 0.0) {
        flat[i] = std::max(flat[i], old_flat);
        fluid_v[FLT][i] = std::max(flat[i] - 1.0, 0.0);
      } else {
        fluid_v[FLT][i] = std::max(fluid_v[FLT][i] - 1.0, 0.0);
        flat[i] = old_flat;
      }
    }
  }

  constexpr inline void evolve() {
    grid.template build_volume<0>();

    auto i_ = Eigen::seqN(frame::i0 + 3, Eigen::fix<frame::iN - 5>);
    dm(i_) = fluid_v[RHO](i_) * grid.dvol[0](i_);
    dtdm(i_) = hydro.dt / dm(i_);
    grid.xe[1](i_) = grid.xe[0](i_) + hydro.dt * umid(i_) / grid.radius;
    upmid(i_) = umid(i_) * pmid(i_);

    for (index_type i = 0; i < 2; ++i) {
      grid.xe[1][frame::i0 + 2 - i] = grid.xe[0][frame::i0 + 2 - i];
      grid.xe[1][frame::iN - 2 + i] = grid.xe[0][frame::iN - 2 + i];
    }

    for (index_type i = frame::i0; i < frame::iN - 1; ++i) {
      grid.dx[1][i] = grid.xe[1][i + 1] - grid.xe[1][i];
      grid.xc[1][i] = grid.xe[1][i] + 0.5 * grid.dx[1][i];
    }

    grid.template build_volume<1>();
    grid.build_amid();

    auto j_ = Eigen::seqN(frame::i0 + 3, Eigen::fix<frame::iN - 6>);

    fluid_v[RHO](j_) *= (grid.dvol[0](j_) / grid.dvol[1](j_));

    uold(j_) = fluid_v[XVL](j_);

    fluid_v[XVL](j_) -= dtdm(j_) * (pmid(j_ + 1) - pmid(j_)) * half *
                        (grid.amid(j_ + 1) + grid.amid(j_));

    fluid_v[ENT](j_) -= dtdm(j_) * (grid.amid(j_ + 1) * upmid(j_ + 1) -
                                    grid.amid(j_) * upmid(j_));

        build_ekin();

    fluid_v[EIN](j_) = fluid_v[ENT](j_) - ekin(j_);
    fluid_v[EIN](j_) =
        fluid_v[EIN](j_).max(smallp / (hydro.gammam1 * fluid_v[RHO](j_)));

  }

  constexpr inline void remap() {
    double fn1 = 0., fn2 = 0.;

    parabolas.prepare(grid.dx[1]);

    parabolas.build(fluid_v, flat, CONSERVATIVE_VARS);

    grid.build_overlap();

    for (index_type f = 0; f < N_FLUIDV; ++f)
      fluxes[f].fill(0.0);

    for (index_type i = frame::j0; i < frame::jM + 1; ++i) {
      if (grid.xndiff[i] >= 0.0) {
        fn1 = 0.5 * grid.xndiff[i] / grid.dx[1][i - 1];
        fn2 = 1.0 - fourthird * fn1;

        fluxes[RHO][i] =
            parabolas[RHO].rc_left(i - 1, fn1, fn2) * grid.xnolap[i];
        // fix this!
        for (index_type f = 1; f < N_FLUIDV; ++f) {
          fluxes[f][i] = parabolas[f].rc_left(i - 1, fn1, fn2) * fluxes[RHO][i];
        }
      } else {
        // tripple check
        fn1 = 0.5 * grid.xndiff[i] / grid.dx[1][i];
        fn2 = 1.0 + fourthird * fn1;

        fluxes[RHO][i] = parabolas[RHO].rc_right(i, -fn1, fn2) * grid.xnolap[i];
        // fix this!
        for (index_type f = 1; f < N_FLUIDV; ++f) {
          fluxes[f][i] = parabolas[f].rc_right(i, -fn1, fn2) * fluxes[RHO][i];
        }
      }
    }

    for (index_type i = frame::j0 - 1; i < frame::jM + 1; ++i) {
      dm[i] = fluid_v[RHO][i] * grid.dvol[1][i];
      dm0[i] = (dm[i] + fluxes[RHO][i] - fluxes[RHO][i + 1]);

      fluid_v[RHO][i] = dm0[i] / grid.dvol[0][i];
      fluid_v[RHO][i] = std::max(smallr, fluid_v[RHO][i]);

      dm0[i] = 1. / (fluid_v[RHO][i] * grid.dvol[0][i]);

      for (index_type f = 1; f < N_FLUIDV; ++f) {
        fluid_v[f][i] =
            (fluid_v[f][i] * dm[i] + fluxes[f][i] - fluxes[f][i + 1]) * dm0[i];
      }
    }

    // supersonic flow
    build_ekin();
    for (index_type i = frame::j0; i < frame::jM; ++i) {
      if (ekin[i] / fluid_v[EIN][i] < 100.0)
        fluid_v[EIN][i] = fluid_v[ENT][i] - ekin[i];
      fluid_v[PRS][i] = hydro.gammam1 * fluid_v[RHO][i] * fluid_v[EIN][i];
      fluid_v[PRS][i] = std::max(smallp, fluid_v[PRS][i]);
    }

  }  // remap

  constexpr inline void wiggle() {
    if (xwig * flat.sum() != 0.0) {
      //      std::cout << "DO WHIG\n";
      grid.store_coords();

      //      for (size_t i = NZ_GHOST + 1; i < nx + NZ_GHOST; ++i) {
      for (index_type i = frame::j0 + 1; i < frame::jM; ++i) {
        if (std::max(flat[i - 1], flat[i]) > 0.0) {
          grid.xe[0][i] += xwig * grid.dx[0][i];
          grid.dx[0][i - 1] = grid.xe[0][i] - grid.xe[0][i - 1];
        }
      }

      grid.template build_volume<0>();
      remap();

      grid.retieve_coords();
    }
  }

  constexpr inline void operator()() {
    bc.apply(grid, fluid_v);

    build_cdx();

    flatten();

    parabolas.prepare(grid.dx[0]);

    parabolas.build(fluid_v, flat, PRIMATIVE_VARS);

    states.prepare(hydro, cdx);

    states.build(parabolas.arr());

    riemann_solve(hydro, states, pmid, umid);

    evolve();

    wiggle();

    remap();
  }
};

}  // namespace ppm
