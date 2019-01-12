#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <utility>

#include "global_def.hpp"
#include "physics/hydro.hpp"
#include "ppmlr/ppm_constant.hpp"
#include "ppmlr/state.hpp"

namespace ppm {

// constexpr index_type N_STATEV = 3;
constexpr double NewRapTol = 1.0E-5;
constexpr index_type NewRapItr = 12;

template <typename StateType>
struct RiemannSolve {
  using state_type = StateType;
  using frame = typename state_type::frame;
  using vector_type = typename state_type::vector_type;

  using state_set = StateSet<StateType, N_STATEV>;

  double gamfac1, gamfac2;

  vector_type cl, cr;
  vector_type pli, pri, rli, rri;
  vector_type wl, wr, zl, zr;
  vector_type umidl, umidr, pmold;

  RiemannSolve() = default;
  ~RiemannSolve() = default;

  inline void operator()(const physics::Hydro& hydro,
                         const state_set& states,
                         vector_type& pmid,
                         vector_type& umid) {
    double gamfac1 = hydro.gamma + one;
    double gamfac2 = half * gamfac1 / hydro.gamma;

    auto i_ = Eigen::seqN(frame::i0 + 3, Eigen::fix<frame::iN - 5>);
    //        auto lo = NZ_GHOST - 3, hi = nx + 2 * NZ_GHOST - 2;

    cl(i_) = (hydro.gamma * states[PRS].l(i_) * states[RHO].l(i_)).sqrt();
    pli(i_) = states[PRS].l(i_).inverse();
    rli(i_) = states[RHO].l(i_).inverse();

    cr(i_) = (hydro.gamma * states[PRS].r(i_) * states[RHO].r(i_)).sqrt();
    pri(i_) = states[PRS].r(i_).inverse();
    rri(i_) = states[RHO].r(i_).inverse();
    /*for (index_type i = frame::i0 + 3; i < frame::iN - 2; ++i) {
      cl[i] = std::sqrt(hydro.gamma * states[PRS].l[i] * states[RHO].l[i]);
      pli[i] = one / states[PRS].l[i];
      rli[i] = one / states[RHO].l[i];
    }*/

    /*for (index_type i = frame::i0 + 3; i < frame::iN - 2; ++i) {
      cr[i] = std::sqrt(hydro.gamma * states[PRS].r[i] * states[RHO].r[i]);
      pri[i] = one / states[PRS].r[i];
      rri[i] = one / states[RHO].r[i];
    }*/
    pmid(i_) = states[PRS].r(i_) - states[PRS].l(i_) -
               cr(i_) * (states[XVL].r(i_) - states[XVL].l(i_));
    pmid(i_) = states[PRS].l(i_) + pmid(i_) * cl(i_) / (cl(i_) + cr(i_));
    pmid(i_) = pmid(i_).max(small_riemann);
    // pmid(i_) = std::max(small_riemann, pmid[i]);
    /*
    for (index_type i = frame::i0 + 3; i < frame::iN - 2; ++i) {
      pmid[i] = states[PRS].r[i] - states[PRS].l[i] -
                cr[i] * (states[XVL].r[i] - states[XVL].l[i]);
      pmid[i] = states[PRS].l[i] + pmid[i] * cl[i] / (cl[i] + cr[i]);
      pmid[i] = std::max(small_riemann, pmid[i]);
    }*/

    for (index_type i = frame::i0 + 3; i < frame::iN - 2; ++i) {
      for (index_type n = 0; n < NewRapItr; ++n) {
        pmold[i] = pmid[i];
        wl[i] = one + gamfac2 * (pmid[i] - states[PRS].l[i]) * pli[i];
        wr[i] = one + gamfac2 * (pmid[i] - states[PRS].r[i]) * pri[i];

        wl[i] = cl[i] * std::sqrt(wl[i]);
        wr[i] = cr[i] * std::sqrt(wr[i]);

        zl[i] = 4.0 * rli[i] * wl[i] * wl[i];
        zr[i] = 4.0 * rri[i] * wr[i] * wr[i];

        zl[i] =
            -zl[i] * wl[i] / (zl[i] - gamfac1 * (pmid[i] - states[PRS].l[i]));
        zr[i] =
            zr[i] * wr[i] / (zr[i] - gamfac1 * (pmid[i] - states[PRS].r[i]));

        umidl[i] = states[XVL].l[i] - (pmid[i] - states[PRS].l[i]) / wl[i];
        umidr[i] = states[XVL].r[i] + (pmid[i] - states[PRS].r[i]) / wr[i];

        pmid[i] += (umidr[i] - umidl[i]) * (zl[i] * zr[i]) / (zr[i] - zl[i]);
        pmid[i] = std::max(small_riemann, pmid[i]);
        if (std::abs(pmid[i] - pmold[i]) / pmid[i] < NewRapTol)
          break;
      }
    }
    umidl(i_) = states[XVL].l(i_) - (pmid(i_) - states[PRS].l(i_)) / wl(i_);
    umidr(i_) = states[XVL].r(i_) + (pmid(i_) - states[PRS].r(i_)) / wr(i_);
    umid(i_) = half * (umidl(i_) + umidr(i_));

    /*    for (index_type i = frame::i0 + 3; i < frame::iN - 2; ++i) {
          umidl[i] = states[XVL].l[i] - (pmid[i] - states[PRS].l[i]) / wl[i];
          umidr[i] = states[XVL].r[i] + (pmid[i] - states[PRS].r[i]) / wr[i];
          umid[i] = half * (umidl[i] + umidr[i]);
        }*/
  }
};

}  // namespace ppm
