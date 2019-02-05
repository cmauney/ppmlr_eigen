#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>

#include <Eigen/Dense>

#include "global_def.hpp"
#include "physics/timestep.hpp"
#include "physics/eos.hpp"
#include "ppmlr/boundary_conditions.hpp"
#include "ppmlr/parabola.hpp"
#include "ppmlr/ppm_constant.hpp"
#include "ppmlr/ppm_def.hpp"
#include "ppmlr/ppm_grid.hpp"
#include "ppmlr/riemann_solve.hpp"
#include "ppmlr/state.hpp"

namespace ppm {

template <typename GeoType, index_type N>
struct PPMLR {
  using frame = VectorFrame<N>;
  using fluid_v = PPMVarSet<N>;

  using grid_type = PPMGrid<GeoType, frame>;
  using bc_type = BoundaryCondition<grid_type>;
  using parabola_basis_type = ParabolaBasis<frame>;
  using parabola_type = Parabola<parabola_basis_type>;
  using state_domain_type = StateDomain<frame>;
  using state_type = State<parabola_type, state_domain_type>;

 // template<index_type NP>
//  using parabola_set = ParabolaSet<parabola_type, NP>;
  using state_set = std::array<state_type, NQV>;
  using riemann = RiemannSolve<state_type>;

  using vector_type = VectorType<frame::iN>; 

//  const physics::Hydro& hydro;  // reference to 'globals'
  const physics::Timestep& ts;
  const physics::IdealGas& eos;

  // max svel in sweep
  double max_cdx;
  // coordinate vars
  grid_type grid;
  // fluid vars
//  FluidVars<vector_type> fluid_v;
  fluid_v v;
  PPMSub<N, NQV> Q;
  PPMSub<N, NUV> U;

  // boundary do-er
  bc_type bc;


  // flatten stuff
  //  flat_type flat;
  // parabolas
  parabola_basis_type             parabola_basis;
  std::array<parabola_type, NQV>  parabolas_prim;
  std::array<parabola_type, NUV>  parabolas_cons;
  // states
  state_domain_type               state_domain;
  std::array<state_type, NQV>     states_prim;
//  state_set states;
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
      : ts(_ts), eos(_eos), max_cdx(zero), 
      Q{v.col(RHO), 
        v.col(PRS), 
        v.col(XVL)
      }, 
      U{v.col(RHO), 
        v.col(XVL), 
        v.col(YVL), 
        v.col(ZVL), 
        v.col(EIN), 
        v.col(ENT)}
  {}

  ~PPMLR() = default;

    template <class V>
    inline void put_vars( const V& xe, const V& xc, const V& dx, const V& rho, const V& prs, const V& xvl, const V& yvl, const V& zvl, const V& f)
    {
    auto i_lhs = Eigen::seqN(frame::j0, Eigen::fix<N>);

    grid.xe0(i_lhs) = xe;
    grid.xc0(i_lhs) = xc;
    grid.dx0(i_lhs) = dx;

    v.col(RHO)(i_lhs) = rho;
    v.col(PRS)(i_lhs) = prs;
    v.col(XVL)(i_lhs) = xvl;
    v.col(YVL)(i_lhs) = yvl;
    v.col(ZVL)(i_lhs) = zvl;
    v.col(FLT)(i_lhs) = f;
    
    build_ekin();
    
    v.col(ENT)(i_lhs) =
        v.col(PRS)(i_lhs) / (v.col(RHO)(i_lhs) * eos.gammam1) +
        ekin(i_lhs);
  }

  template <class V>
  inline void get_vars(V& rho,
                       V& prs,
                       V& xvel,
                       V& yvel,
                       V& zvel,
                       V& f) {

    auto i_rhs = Eigen::seqN(frame::j0, Eigen::fix<N>);
    
    rho   = v.col(RHO)(i_rhs);
    prs   = v.col(PRS)(i_rhs);
    xvel  = v.col(XVL)(i_rhs);
    yvel  = v.col(YVL)(i_rhs);
    zvel  = v.col(ZVL)(i_rhs);
    f     = v.col(FLT)(i_rhs);
  }

  constexpr inline void build_ekin() {
    ekin = half * (v.col(XVL).square() + v.col(YVL).square() +
                   v.col(ZVL).square());
  }

  constexpr inline void build_cdx() {
    auto i_ = Eigen::seqN(frame::i0 + 2, Eigen::fix<frame::iN - 4>);

    cdx(i_) = (eos.gamma * v.col(PRS)(i_) / v.col(RHO)(i_)).sqrt() /
              (grid.dx0(i_) * grid.radius);

    max_cdx = cdx(i_).maxCoeff();
  }

  inline void flatten() {
    double delp1, delp2, zsten;
    double shock, old_flat;

    for (index_type i = frame::j0 - 4; i < frame::jM + 4; ++i) {
      delp1 = v.col(PRS)[i + 1] - v.col(PRS)[i - 1];
      delp2 = v.col(PRS)[i + 2] - v.col(PRS)[i - 2];
      if (std::abs(delp2) < small)
        delp2 = small;
      shock =
          std::abs(delp1) / std::min(v.col(PRS)[i + 1], v.col(PRS)[i - 1]) -
          epsilon;
      shock = std::max(0.0, shock);
      if (shock > 0.0)
        shock = 1.0;
      if (v.col(XVL)[i - 1] < v.col(XVL)[i + 1])
        shock = 0.0;

      z[i] = shock * std::max(0.0, (delp1 / delp2 - omega1) * omega2);
    }

    z[1] = z[2];
    z[frame::iN - 2] = z[frame::iN - 3];

    flat = zero;
    for (index_type i = frame::j0 - 4; i < frame::jM + 4; ++i) {
      zsten = std::max(z[i - 1], std::max(z[i], z[i + 1]));
      flat[i] = std::max(0.0, std::min(0.5, zsten));
    }

    for (index_type i = frame::j0 - 3; i < frame::jM + 3; ++i) {
      old_flat = v.col(FLT)[i] - static_cast<int>(v.col(FLT)[i]);
      if (flat[i] > 0.0) {
        flat[i] = std::max(flat[i], old_flat);
        v.col(FLT)[i] = std::max(flat[i] - 1.0, 0.0);
      } else {
        v.col(FLT)[i] = std::max(v.col(FLT)[i] - 1.0, 0.0);
        flat[i] = old_flat;
      }
    }
  }

  constexpr inline void evolve() {
    grid.template build_volume<0>();

    auto i_ = Eigen::seqN(frame::i0 + 3, Eigen::fix<frame::iN - 5>);
    dm(i_) = v.col(RHO)(i_) * grid.dvol0(i_);
    dtdm(i_) = ts.dt / dm(i_);
    grid.xen(i_) = grid.xe0(i_) + ts.dt * umid(i_) / grid.radius;
    upmid(i_) = umid(i_) * pmid(i_);

    for (index_type i = 0; i < 2; ++i) {
      grid.xen[frame::i0 + 2 - i] = grid.xe0[frame::i0 + 2 - i];
      grid.xen[frame::iN - 2 + i] = grid.xe0[frame::iN - 2 + i];
    }

    for (index_type i = frame::i0; i < frame::iN - 1; ++i) {
      grid.dxn[i] = grid.xen[i + 1] - grid.xen[i];
      grid.xcn[i] = grid.xen[i] + 0.5 * grid.dxn[i];
    }

    grid.template build_volume<1>();
    grid.build_amid();

    auto j_ = Eigen::seqN(frame::i0 + 3, Eigen::fix<frame::iN - 6>);

    v.col(RHO)(j_) *= (grid.dvol0(j_) / grid.dvoln(j_));

    uold(j_) = v.col(XVL)(j_);

    v.col(XVL)(j_) -= dtdm(j_) * (pmid(j_ + 1) - pmid(j_)) * half *
                        (grid.amid(j_ + 1) + grid.amid(j_));

    v.col(ENT)(j_) -= dtdm(j_) * (grid.amid(j_ + 1) * upmid(j_ + 1) -
                                    grid.amid(j_) * upmid(j_));

        build_ekin();

    v.col(EIN)(j_) = v.col(ENT)(j_) - ekin(j_);
    v.col(EIN)(j_) =
        v.col(EIN)(j_).max(smallp / (eos.gammam1 * v.col(RHO)(j_)));

  }

  constexpr inline void remap() {
    double fn1 = 0., fn2 = 0.;

//    parabolas.prepare(grid.dx[1]);
    parabola_basis.prepare(grid.dxn); 

    for(index_type i = 0; i < NUV; ++i)
      parabolas_cons[i].construct(U[i], flat, parabola_basis.basis);

    grid.build_overlap();

    for (index_type f = 0; f < NUV; ++f)
      fluxes[f].fill(0.0);

    for (index_type i = frame::j0; i < frame::jM + 1; ++i) {
      if (grid.xndiff[i] >= 0.0) {
        fn1 = 0.5 * grid.xndiff[i] / grid.dxn[i - 1];
        fn2 = 1.0 - fourthird * fn1;

        fluxes[URHO][i] =
            parabolas_cons[URHO].rc_left(i - 1, fn1, fn2) * grid.xnolap[i];
        // fix this!
        for (index_type f = 1; f < NUV; ++f) {
          fluxes[f][i] = parabolas_cons[f].rc_left(i - 1, fn1, fn2) * fluxes[URHO][i];
        }
      } else {
        // tripple check
        fn1 = 0.5 * grid.xndiff[i] / grid.dxn[i];
        fn2 = 1.0 + fourthird * fn1;

        fluxes[URHO][i] = parabolas_cons[URHO].rc_right(i, -fn1, fn2) * grid.xnolap[i];
        // fix this!
        for (index_type f = 1; f < NUV; ++f) {
          fluxes[f][i] = parabolas_cons[f].rc_right(i, -fn1, fn2) * fluxes[URHO][i];
        }
      }
    }

    for (index_type i = frame::j0 - 1; i < frame::jM + 1; ++i) {
      dm[i] = U[URHO][i] * grid.dvoln[i];
      dm0[i] = (dm[i] + fluxes[URHO][i] - fluxes[URHO][i + 1]);

      U[URHO][i] = dm0[i] / grid.dvol0[i];
      U[URHO][i] = std::max(smallr, U[URHO][i]);

      dm0[i] = 1. / (U[URHO][i] * grid.dvol0[i]);

      for (index_type f = 1; f < NUV; ++f) {
        U[f][i] =
            (U[f][i] * dm[i] + fluxes[f][i] - fluxes[f][i + 1]) * dm0[i];
      }
    }

    // supersonic flow
    build_ekin();
    for (index_type i = frame::j0; i < frame::jM; ++i) {
      if (ekin[i] / U[EIN][i] < 100.0)
        U[UEIN][i] = U[ENT][i] - ekin[i];
      v.col(PRS)[i] = eos.gammam1 * U[URHO][i] * U[UEIN][i];
      v.col(PRS)[i] = std::max(smallp, v.col(PRS)[i]);
    }

  }  // remap

  constexpr inline void wiggle() {
    if (xwig * flat.sum() != 0.0) {
      //      std::cout << "DO WHIG\n";
      grid.store_coords();

      //      for (size_t i = NZ_GHOST + 1; i < nx + NZ_GHOST; ++i) {
      for (index_type i = frame::j0 + 1; i < frame::jM; ++i) {
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

  constexpr inline void operator()() {
    bc.apply(grid, v);

    build_cdx();

    flatten();

    parabola_basis.prepare(grid.dx0);
    state_domain.build(ts.hdt, cdx);

    for(index_type i = 0; i < NQV; ++i){
      parabolas_prim[i].construct(Q[i], flat, parabola_basis.basis);
      states_prim[i].build_state(parabolas_prim[i], state_domain);
    }


    riemann_solve(eos, states_prim, pmid, umid);

    evolve();

    wiggle();

    remap();
  }
};

}  // namespace ppm
