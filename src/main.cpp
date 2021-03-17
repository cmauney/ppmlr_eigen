#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>

#include <fmt/core.h>
#include <string_view>

#include "physics/hydro.hpp"
#include "ppmlr/ppm_grid.hpp"
#include "ppmlr/ppmlr.hpp"
#include "ppmlr/riemann_solve.hpp"
#include "simulation.hpp"

//#define _GNU_SOURCE
//#include <fenv.h>

constexpr index_type NCYC = 100;

constexpr double pr = 0.1;
constexpr double dr = 0.125;
constexpr double pl = 1.0;
constexpr double dl = 1.0;

template<typename V>
constexpr inline auto
etot(index_type n, const V& p, const V& r, const std::array<V, 3>& vs, V& e)
{
  double gm1 = (simulation::my_gamma - 1.0);
  e = p / (r * gm1) + half * (vs[0].square() + vs[1].square() + vs[2].square());
  //  for (index_type i = 0; i < n; ++i) {
  //    e[i] =
  //        p[i] / (r[i] * gm1) +
  //        0.5 * (vs[0][i] * vs[0][i] + vs[1][i] * vs[1][i] + vs[2][i] *
  //        vs[2][i]);
  //  }
}

template<template<index_type> class V, index_type N>
auto
gridder(double xmin, double xmax)
{
  V<N> xe, xc, dx;

  auto ddx = (xmax - xmin) / static_cast<double>(simulation::NX);

  for (index_type i = 0; i < simulation::NX; ++i) {
    dx[i] = ddx;
    xe[i] = xmin + static_cast<double>(i) * dx[i];
    xc[i] = xe[i] + 0.5 * dx[i];
  }

  return std::make_tuple(xe, xc, dx);
}

template<template<index_type> class V, index_type N>
auto
initial_conditions(const V<N>& xc)
{
  V<N> rho, prs, ene, fla;
  std::array<V<N>, 3> vels;

  for (index_type i = 0; i < N; ++i) {
    fla[i] = 0.0;
    if (xc[i] <= 0.5) {
      rho[i] = dl;
      prs[i] = pl;
    } else {
      rho[i] = dr;
      prs[i] = pr;
    }
  }

  for (auto& v : vels) {
    for (index_type i = 0; i < N; ++i) {
      v[i] = 0.0;
    }
  }

  etot(N, prs, rho, vels, ene);

  return std::make_tuple(rho, prs, vels, ene, fla);
}

template<typename... Vs>
void
dump(index_type N, double t, Vs&&... vs)
{
  std::ofstream of("run00.dat", std::ofstream::out | std::ofstream::app);

  std::string s;
  auto commd = [](double a) { return std::to_string(a) + ", "; };
  of << t << std::endl;
  for (index_type i = 0; i < N; ++i) {
    s = "";
    s = (commd(vs[i]) + ...);
    of << i << " " << s << std::endl;
  }
  of << "\n\n";
  of.close();
}

struct sweep_x
{
  ppm::PPMLR<simulation::XGEO, simulation::NX> ppmlr_drv;

  sweep_x(const physics::Timestep& ts, const physics::IdealGas& eos)
    : ppmlr_drv(ts, eos)
  {}
};

int
main()
{
  //  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  //  physics::Hydro hydro(my_gamma);
  physics::IdealGas eos(simulation::my_gamma);
  physics::Timestep ts;
  sweep_x swpx(ts, eos);

  auto [xe, xc, dx] = gridder<VectorType, simulation::NX>(0.0, 1.0);
  auto [r, p, u, e, f] = initial_conditions<VectorType, simulation::NX>(xc);

  std::chrono::duration<double> rt_s(0);

  ts.first_dt(eos, p, r, u[0], dx);
  double time = 0.0;

  index_type ncycle = 0;
  for (index_type it = 0; it < NCYC; ++it) {
    dump(simulation::NX, time, xe, r, p, u[0], e, f);

    auto start_hydro = std::chrono::high_resolution_clock::now();

    time += ts.dt;
    ++ncycle;

    etot(simulation::NX, p, r, u, e);

    auto end_hydro = std::chrono::high_resolution_clock::now();
    rt_s += (end_hydro - start_hydro);

    swpx.ppmlr_drv.put_vars(xe, xc, dx, r, p, u[0], u[1], u[2], f);

    start_hydro = std::chrono::high_resolution_clock::now();
    //    sweepx();
    swpx.ppmlr_drv();

    end_hydro = std::chrono::high_resolution_clock::now();
    rt_s += (end_hydro - start_hydro);

    swpx.ppmlr_drv.get_vars(r, p, u[0], u[1], u[2], f);

    start_hydro = std::chrono::high_resolution_clock::now();

    ts.dt_courant(swpx.ppmlr_drv.max_cdx, u[0], dx);

    end_hydro = std::chrono::high_resolution_clock::now();
    rt_s += (end_hydro - start_hydro);
  }

  std::cout << time << std::endl;

  //    std::cout << "ttf = " << hydro_duration.count() << " us\n";
  std::cout << "ttf = " << rt_s.count() * 1.0E-6 << " us\n";

  return 0;
}
