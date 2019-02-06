#pragma once

#include <array>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>

#include "global_def.hpp"
#include "ppm_constant.hpp"
#include "ppm_def.hpp"
#include "ppm_grid.hpp"

namespace ppm {
template<index_type N>
struct ParabolaBasis
{
  //  using frame = Frame;
  using vector_type = __PPMVector<N>;
  using arg_vector_type = __PPMVecRef<N>;
  using matrix_type = __PPMMatrix<N, 5>;
  using arg_matrix_type = __PPMMatRef<N, 5>;
  //  using arg_vector_type = ArgConstType<vector_type>;

  // using basis_type = fixed_storage_type<bp::iN * 5>;
  //  using matrix_type = MatrixType<frame::iN, 5>;
  //  using arg_matrix_type = ArgConstType<matrix_type>;

  // internal variables
  vector_type a, b, c, d, ai, bi, ci;

  // basis:
  //  std::array<double, bp::iN * 5> basis;
  matrix_type basis;

  // spans/views
  ParabolaBasis() = default;
  ~ParabolaBasis() = default;

  constexpr inline void clear()
  {
    a = zero;
    b = zero;
    c = zero;
    d = zero;
    ai = zero;
    bi = zero;
    ci = zero;
    basis = zero;
  }

  // constexpr inline void prepare(const vector_type& dx) {
  constexpr inline void prepare(const vector_type& dx)
  {
    using Eigen::fix;
    using Eigen::seqN;

    auto i_ = seqN(0, fix<__PPM_N<N> - 1>);
    a(i_) = dx(i_) + dx(i_ + 1);
    b(i_) = a(i_) + dx(i_);
    c(i_) = a(i_) + dx(i_ + 1);
    ai(i_) = 1. / a(i_);
    bi(i_) = 1. / b(i_);
    ci(i_) = 1. / c(i_);

    auto j_ = seqN(1, fix<__PPM_N<N> - 3>);
    d(j_) = 1. / (a(j_ - 1) + a(j_ + 1));
    basis(j_, 0) =
      dx(j_) * ai(j_) + 2.0 * dx(j_ + 1) * dx(j_) * d(j_) * ai(j_) *
                          (a(j_ - 1) * bi(j_) - a(j_ + 1) * ci(j_));

    basis(j_, 1) = -d(j_) * dx(j_) * a(j_ - 1) * bi(j_);
    basis(j_, 2) = d(j_) * dx(j_ + 1) * a(j_ + 1) * ci(j_);

    d(j_) = dx(j_) / (a(j_ - 1) + dx(j_ + 1));
    basis(j_, 3) = d(j_) * b(j_ - 1) * ai(j_);
    basis(j_, 4) = d(j_) * c(j_) * ai(j_ - 1);
    // check this range - upd: should be fine, one extra term worth not
    // restarting loop;
  }
};
template<index_type N>
struct Parabola
{
  using vector_type = __PPMVector<N>;
  using arg_vector_type = __PPMVecRef<N>;
  using matrix_type = __PPMMatrix<N, 5>;
  using arg_matrix_type = __PPMMatRef<N, 5>;

  // parabola coeffs
  vector_type da, a6, al, ar;

  // internals
  vector_type diff_an, del_a;
  vector_type diff_amin;
  vector_type tmp1, tmp2, tmp3, onemfl;

  Parabola() = default;
  ~Parabola() = default;

  constexpr inline void clear()
  {
    diff_an = zero;
    del_a = zero;

    tmp1 = zero;
    tmp2 = zero;
    tmp3 = zero;

    da = zero;
    a6 = zero;
    al = zero;
    ar = zero;
  }

  constexpr inline void construct(const vector_type& a,
                                  const vector_type& flat,
                                  const matrix_type& pbase)
  {
    // clear();

    using Eigen::fix;
    using Eigen::seqN;

    //    constexpr auto i0 = frame::i0;
    //    constexpr auto iN = frame::iN;
    //    constexpr auto nN = iN - i0;

    const auto Nm1 = fix<__PPM_N<N> - 1>;
    const auto Nm2 = fix<__PPM_N<N> - 2>;
    const auto Nm3 = fix<__PPM_N<N> - 3>;
    const auto Nm4 = fix<__PPM_N<N> - 4>;

    auto i_ = seqN(0, Nm1);
    diff_an(i_) = a(i_ + 1) - a(i_);

    auto j_ = seqN(1, Nm2);
    del_a(j_) = pbase(j_, 3) * diff_an(j_) + pbase(j_, 4) * diff_an(j_ - 1);

    diff_amin(j_) = diff_an(j_ - 1).min(diff_an(j_));

    del_a(j_) *= (del_a(j_)).abs().min(two * diff_amin(j_)).sign();
    del_a(j_) = (diff_an(j_ - 1) * diff_an(j_)).select(zero, del_a(j_));

    auto k_ = seqN(1, Nm3);
    ar(k_) = a(k_) + pbase(k_, 0) * diff_an(k_) + pbase(k_, 1) * del_a(k_ + 1) +
             pbase(k_, 2) * del_a(k_);

    al(k_ + 1) = ar(k_);

    auto l_ = seqN(2, Nm4);
    onemfl(l_) = 1.0 - flat(l_);

    ar(l_) = flat(l_) * a(l_) + onemfl(l_) * ar(l_);
    al(l_) = flat(l_) * a(l_) + onemfl(l_) * al(l_);

    da(l_) = ar(l_) - al(l_);
    a6(l_) = (6.0 * (a(l_) - half * (al(l_) + ar(l_))));
    tmp1(l_) = (ar(l_) - a(l_)) * (a(l_) - al(l_));
    tmp2(l_) = da(l_).square();
    tmp3(l_) = da(l_) * a6(l_);

    for (index_type i = 2; i < Nm2; ++i) {
      bool redo = false;
      if (tmp1[i] <= zero) {
        ar[i] = a[i];
        al[i] = a[i];
        redo = true;
      }
      if (tmp2[i] < +tmp3[i]) {
        al[i] = three * a[i] - two * ar[i];
        redo = true;
      }
      if (tmp2[i] < -tmp3[i]) {
        ar[i] = three * a[i] - two * al[i];
        redo = true;
      }
      if (redo) {
        da[i] = ar[i] - al[i];
        a6[i] = 6.0 * (a[i] - half * (al[i] + ar[i]));
      }
    }
  }

  constexpr inline double rc_left(const index_type i,
                                  const double f1,
                                  const double f2) const
  {
    return ar[i] - f1 * (da[i] - f2 * a6[i]);
  }

  constexpr inline double rc_right(const index_type i,
                                   const double f1,
                                   const double f2) const
  {
    return al[i] + f1 * (da[i] + f2 * a6[i]);
  }

  template<index_type IDX0, index_type IDXN>
  constexpr inline decltype(auto) rc_leftV(const vector_type& f1,
                                           const vector_type& f2) const
  {
    auto i_ = Eigen::seqN(IDX0, Eigen::fix<IDXN>);
    return ar(i_) - f1(i_) * (da(i_) - f2(i_) * a6(i_));
  }

  template<index_type IDX0, index_type IDXN>
  constexpr inline decltype(auto) rc_rightV(const vector_type& f1,
                                            const vector_type& f2) const
  {
    auto i_ = Eigen::seqN(IDX0, Eigen::fix<IDXN>);
    return al(i_) + f1(i_) * (da(i_) + f2(i_) * a6(i_));
  }
};

} // namespace ppm
