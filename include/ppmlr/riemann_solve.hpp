#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <utility>

#include "global_def.hpp"
#include "physics/eos.hpp"
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

  using state_set = std::array<StateType, NQV>;

  double gamfac1, gamfac2;

  vector_type cl, cr;
  vector_type pli, pri, rli, rri;
  vector_type wl, wr, zl, zr;
  vector_type umidl, umidr, pmold;

  RiemannSolve() = default;
  ~RiemannSolve() = default;

  inline void operator()(const physics::IdealGas& eos,
                         const state_set& states,
                         vector_type& pmid,
                         vector_type& umid) {
    double gamfac1 = eos.gamma + one;
    double gamfac2 = half * gamfac1 / eos.gamma;

    auto i_ = Eigen::seqN(frame::i0 + 3, Eigen::fix<frame::iN - 5>);

    cl(i_) = (eos.gamma * states[QPRS].l(i_) * states[QRHO].l(i_)).sqrt();
    pli(i_) = states[QPRS].l(i_).inverse();
    rli(i_) = states[QRHO].l(i_).inverse();

    cr(i_) = (eos.gamma * states[QPRS].r(i_) * states[QRHO].r(i_)).sqrt();
    pri(i_) = states[QPRS].r(i_).inverse();
    rri(i_) = states[QRHO].r(i_).inverse();
    pmid(i_) = states[QPRS].r(i_) - states[QPRS].l(i_) -
               cr(i_) * (states[QXVL].r(i_) - states[QXVL].l(i_));
    pmid(i_) = states[QPRS].l(i_) + pmid(i_) * cl(i_) / (cl(i_) + cr(i_));
    pmid(i_) = pmid(i_).max(small_riemann);

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

        umidl[i] = states[QXVL].l[i] - (pmid[i] - states[PRS].l[i]) / wl[i];
        umidr[i] = states[QXVL].r[i] + (pmid[i] - states[PRS].r[i]) / wr[i];

        pmid[i] += (umidr[i] - umidl[i]) * (zl[i] * zr[i]) / (zr[i] - zl[i]);
        pmid[i] = std::max(small_riemann, pmid[i]);
        if (std::abs(pmid[i] - pmold[i]) / pmid[i] < NewRapTol)
          break;
      }
    }
    umidl(i_) = states[QXVL].l(i_) - (pmid(i_) - states[PRS].l(i_)) / wl(i_);
    umidr(i_) = states[QXVL].r(i_) + (pmid(i_) - states[PRS].r(i_)) / wr(i_);
    umid(i_) = half * (umidl(i_) + umidr(i_));

  }
};

}  // namespace ppm
