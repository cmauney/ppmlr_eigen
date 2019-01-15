#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <utility>

#include "global_def.hpp"
#include "parabola.hpp"
#include "physics/hydro.hpp"
#include "ppm_constant.hpp"
#include "ppm_def.hpp"
#include "ppm_grid.hpp"

namespace ppm {
template <typename Frame>
struct StateDomain {
  using frame = Frame;
  using vector_type = VectorType<frame::iN>;

  vector_type Cdtdx, fCdtdx;

  StateDomain() = default;
  ~StateDomain() = default;

  constexpr inline void build(const double hdt,
                              const vector_type& cdx) {
    auto i_ = Eigen::seqN(frame::i0 + 2, Eigen::fix<frame::iN - 4>);
    Cdtdx(i_) = cdx(i_) * hdt;
    fCdtdx(i_) = 1.0 - fourthird * Cdtdx(i_);
    /*
    for (index_type i = frame::i0 + 2; i < frame::iN - 2; ++i) {
      Cdtdx[i] = cdx[i] * hydro.hdt;
      fCdtdx[i] = 1.0 - fourthird * Cdtdx[i];
    }
    */
  }
};

template <typename ReconstructType, typename DomainType>
struct State {
  using frame = typename ReconstructType::frame;
  using vector_type = typename ReconstructType::vector_type;

  using domain_type = DomainType;
  using reconstruct_type = ReconstructType;

  vector_type l, r;
  State() = default;

  constexpr inline void build_state(const reconstruct_type& reco,
                                    const domain_type& domain) {
    for (index_type i = frame::i0 + 2; i < frame::iN - 2; ++i) {
      l[i + 1] = reco.rc_left(i, domain.Cdtdx[i], domain.fCdtdx[i]);
      r[i] = reco.rc_right(i, domain.Cdtdx[i], domain.fCdtdx[i]);
    }
  }
};

template <typename StateType>
struct StateSet {
  using state_type = StateType;
  using vector_type = typename state_type::vector_type;
  using domain_type = typename state_type::domain_type;
  using reconstruct_type = typename state_type::reconstruct_type;

  domain_type domain;

  std::array<state_type, 3> states;

//  state_type rho_s, prs_s, xvl_s;

  StateSet() = default;

  constexpr inline void prepare(const double hdt,
                                const vector_type& cdx) {
    domain.build(hdt, cdx);
  }

  template <typename RCs>
  constexpr inline void build(const RCs& recos) {
    for (index_type i = 0; i < 3; ++i) {
      states[i].build_state(recos[i], domain);
    }
  }

  constexpr inline const state_type& operator[](const index_type i) const {
    return states[i];
  }
};

}  // namespace ppm
