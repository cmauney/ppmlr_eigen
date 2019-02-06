#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>

#include <Eigen/Dense>

#include "global_def.hpp"
#include "physics/eos.hpp"
#include "physics/timestep.hpp"
#include "ppmlr/boundary_conditions.hpp"
#include "ppmlr/parabola.hpp"
#include "ppmlr/ppm_constant.hpp"
#include "ppmlr/ppm_def.hpp"
#include "ppmlr/ppm_grid.hpp"
#include "ppmlr/riemann_solve.hpp"
#include "ppmlr/state.hpp"

namespace ppm {

template<typename GeoType, index_type N>
struct PPMLR
{
  using grid_type = PPMGrid<GeoType, N>;
  using bc_type = BoundaryCondition<grid_type, N>;
  using parabola_basis_type = ParabolaBasis<N>;
  using parabola_type = Parabola<N>;
  using state_domain_type = StateDomain<N>;
  using state_type = State<parabola_type, N>;

  using state_set = std::array<state_type, NQV>;
  using riemann = RiemannSolve<state_type, N>;

  using vector_type = __PPMVector<N>;
  using arg_vector_type = __PPMVecRef<N>;

  const physics::Timestep& ts;
  const physics::IdealGas& eos;

  // max svel in sweep
  double max_cdx;
  // coordinate vars
  grid_type grid;
  // fluid vars
  __PPMFluidVars<N> v;
  __PPMFluidRef<N, NQV> Q;
  __PPMFluidRef<N, NUV> U;

  // boundary do-er
  bc_type bc;

  // parabolas
  parabola_basis_type parabola_basis;
  std::array<parabola_type, NQV> parabolas_prim;
  std::array<parabola_type, NUV> parabolas_cons;
  // states
  state_domain_type state_domain;
  std::array<state_type, NQV> states_prim;
  // riemann-er
  riemann riemann_solve;

  vector_type flat, z;
  // kinetic E store, cspeed
  vector_type ekin, cdx;
  // evolve stores
  vector_type dm, dm0, dtdm, pmid, umid, upmid, uold;
  // remap stores
  std::array<vector_type, NUV> fluxes;

  PPMLR(const physics::Timestep& _ts, const physics::IdealGas& _eos)
    : ts(_ts)
    , eos(_eos)
    , max_cdx(zero)
    , Q{ v[RHO], v[PRS], v[XVL] }
    , U{ v[RHO], v[XVL], v[YVL], v[ZVL], v[EIN], v[ENT] }
  {}

  ~PPMLR() = default;

  template<class V>
  constexpr inline void put_vars(const V& xe,
                                 const V& xc,
                                 const V& dx,
                                 const V& rho,
                                 const V& prs,
                                 const V& xvl,
                                 const V& yvl,
                                 const V& zvl,
                                 const V& f)
  {
    // auto i_lhs = Eigen::seqN(frame::j0, Eigen::fix<N>);
    auto i_lhs = Eigen::seqN(NZ_GHOST, Eigen::fix<N>);

    grid.set_g0(xe, xc, dx);

    v[RHO](i_lhs) = rho;
    v[PRS](i_lhs) = prs;
    v[XVL](i_lhs) = xvl;
    v[YVL](i_lhs) = yvl;
    v[ZVL](i_lhs) = zvl;
    v[FLT](i_lhs) = f;

    build_ekin();

    v[ENT](i_lhs) = v[PRS](i_lhs) / (v[RHO](i_lhs) * eos.gammam1) + ekin(i_lhs);
  }

  template<class V>
  constexpr inline void get_vars(V& rho,
                                 V& prs,
                                 V& xvel,
                                 V& yvel,
                                 V& zvel,
                                 V& f)
  {

    auto i_rhs = Eigen::seqN(NZ_GHOST, Eigen::fix<N>);

    rho = v[RHO](i_rhs);
    prs = v[PRS](i_rhs);
    xvel = v[XVL](i_rhs);
    yvel = v[YVL](i_rhs);
    zvel = v[ZVL](i_rhs);
    f = v[FLT](i_rhs);
  }

  constexpr inline void build_ekin()
  {
    ekin = half * (v[XVL].square() + v[YVL].square() + v[ZVL].square());
  }

  constexpr inline void build_cdx()
  {
    auto i_ = Eigen::seqN(2, Eigen::fix<__PPM_N<N> - 4>);
    cdx(i_) = (eos.gamma * v[PRS](i_) / v[RHO](i_)).sqrt() /
              (grid.dx0(i_) * grid.radius);
    max_cdx = cdx(i_).maxCoeff();
  }

  constexpr inline void flatten()
  {
    double delp1(0), delp2(0), zsten(0);
    double shock(0), old_flat(0);

    for (index_type i = 2; i < __PPM_N<N> - 2; ++i) {
      delp1 = v[PRS][i + 1] - v[PRS][i - 1];
      delp2 = v[PRS][i + 2] - v[PRS][i - 2];
      if (std::abs(delp2) < small)
        delp2 = small;
      shock =
        std::abs(delp1) / std::min(v[PRS][i + 1], v[PRS][i - 1]) - epsilon;
      shock = std::max(0.0, shock);
      if (shock > 0.0)
        shock = 1.0;
      if (v[XVL][i - 1] < v[XVL][i + 1])
        shock = 0.0;

      z[i] = shock * std::max(0.0, (delp1 / delp2 - omega1) * omega2);
    }

    z[1] = z[2];
    z[__PPM_N<N> - 2] = z[__PPM_N<N> - 3];

    flat = zero;
    for (index_type i = 2; i < __PPM_N<N> - 2; ++i) {
      zsten = std::max(z[i - 1], std::max(z[i], z[i + 1]));
      flat[i] = std::max(0.0, std::min(0.5, zsten));
    }

    for (index_type i = 3; i < __PPM_N<N> - 3; ++i) {
      old_flat = v[FLT][i] - static_cast<int>(v[FLT][i]);
      if (flat[i] > 0.0) {
        flat[i] = std::max(flat[i], old_flat);
        v[FLT][i] = std::max(flat[i] - 1.0, 0.0);
      } else {
        v[FLT][i] = std::max(v[FLT][i] - 1.0, 0.0);
        flat[i] = old_flat;
      }
    }
  }

  constexpr inline void evolve()
  {
    grid.template build_volume<0>();

    auto i_ = Eigen::seqN(3, Eigen::fix<__PPM_N<N> - 5>);
    dm(i_) = v[RHO](i_) * grid.dvol0(i_);
    dtdm(i_) = ts.dt / dm(i_);
    grid.xen(i_) = grid.xe0(i_) + ts.dt * umid(i_) / grid.radius;
    upmid(i_) = umid(i_) * pmid(i_);

    for (index_type i = 0; i < 2; ++i) {
      grid.xen[2 - i] = grid.xe0[2 - i];
      grid.xen[__PPM_N<N> - 2 + i] = grid.xe0[__PPM_N<N> - 2 + i];
    }

    auto j_ = Eigen::seqN(0, Eigen::fix<__PPM_N<N> - 1>);
    grid.dxn(j_) = grid.xen(j_ + 1) - grid.xen(j_);
    grid.xcn(j_) = grid.xen(j_) + half * grid.dxn(j_);
    grid.template build_volume<1>();
    grid.build_amid();

    auto k_ = Eigen::seqN(3, Eigen::fix<__PPM_N<N> - 6>);

    v[RHO](k_) *= (grid.dvol0(k_) / grid.dvoln(k_));
    uold(k_) = v[XVL](k_);
    v[XVL](k_) -= dtdm(k_) * (pmid(k_ + 1) - pmid(k_)) * half *
                  (grid.amid(k_ + 1) + grid.amid(k_));
    v[ENT](k_) -= dtdm(k_) * (grid.amid(k_ + 1) * upmid(k_ + 1) -
                              grid.amid(k_) * upmid(k_));
    build_ekin();
    v[EIN](k_) = v[ENT](k_) - ekin(k_);
    v[EIN](k_) = v[EIN](k_).max(smallp / (eos.gammam1 * v[RHO](k_)));
  }

  constexpr inline void remap()
  {
    double fn1 = 0., fn2 = 0.;

    parabola_basis.prepare(grid.dxn);

    for (index_type i = 0; i < NUV; ++i)
      parabolas_cons[i].construct(U[i], flat, parabola_basis.basis);

    grid.build_overlap();

    for (index_type f = 0; f < NUV; ++f)
      fluxes[f].fill(0.0);


    for (index_type i = NZ_GHOST; i < NZ_GHOST + N + 1; ++i) {
      if (grid.xndiff[i] >= 0.0) {
        fn1 = 0.5 * grid.xndiff[i] / grid.dxn[i - 1];
        fn2 = 1.0 - fourthird * fn1;

        fluxes[URHO][i] =
          parabolas_cons[URHO].rc_left(i - 1, fn1, fn2) * grid.xnolap[i];
        // fix this!
        for (index_type f = 1; f < NUV; ++f) {
          fluxes[f][i] =
            parabolas_cons[f].rc_left(i - 1, fn1, fn2) * fluxes[URHO][i];
        }
      } else {
        // tripple check
        fn1 = 0.5 * grid.xndiff[i] / grid.dxn[i];
        fn2 = 1.0 + fourthird * fn1;

        fluxes[URHO][i] =
          parabolas_cons[URHO].rc_right(i, -fn1, fn2) * grid.xnolap[i];
        // fix this!
        for (index_type f = 1; f < NUV; ++f) {
          fluxes[f][i] =
            parabolas_cons[f].rc_right(i, -fn1, fn2) * fluxes[URHO][i];
        }
      }
    }

    for (index_type i = NZ_GHOST - 1; i < NZ_GHOST + N + 1; ++i) {
      dm[i] = U[URHO][i] * grid.dvoln[i];
      dm0[i] = (dm[i] + fluxes[URHO][i] - fluxes[URHO][i + 1]);
      
      U[URHO][i] = dm0[i] / grid.dvol0[i];
      U[URHO][i] = std::max(smallr, U[URHO][i]);
      
      dm0[i] = 1. / (U[URHO][i] * grid.dvol0[i]);
      for (index_type f = 1; f < NUV; ++f) 
        U[f][i] = (U[f][i] * dm[i] + fluxes[f][i] - fluxes[f][i + 1]) * dm0[i];
      
    }

    // supersonic flow
    build_ekin();
    for (index_type i = NZ_GHOST; i < N + NZ_GHOST; ++i) {
      if (ekin[i] / U[EIN][i] < 100.0)
        U[UEIN][i] = U[ENT][i] - ekin[i];
      v[PRS][i] = eos.gammam1 * U[URHO][i] * U[UEIN][i];
      v[PRS][i] = std::max(smallp, v[PRS][i]);
    }

  } // remap

  constexpr inline void wiggle()
  {
    if (xwig * flat.sum() != 0.0) {
      grid.store_coords();

      for (index_type i = NZ_GHOST + 1; i < NZ_GHOST + N; ++i) {
        if (std::max(flat[i - 1], flat[i]) > 0.0) {
          grid.xe0[i] += xwig * grid.dx0[i];
          grid.dx0[i - 1] = grid.xe0[i] - grid.xe0[i - 1];
        }
      }

      grid.template build_volume<0>();
      remap();
      grid.retieve_coords();
    }
  }

  constexpr inline void operator()()
  {
    bc.apply(grid, v);
    build_cdx();
    flatten();

    parabola_basis.prepare(grid.dx0);
    state_domain.build(ts.hdt, cdx);

    for (index_type i = 0; i < NQV; ++i) {
      parabolas_prim[i].construct(Q[i], flat, parabola_basis.basis);
      states_prim[i].build_state(parabolas_prim[i], state_domain);
    }

    riemann_solve(eos, states_prim, pmid, umid);
    evolve();
    wiggle();
    remap();
  }
};

} // namespace ppm
