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
template <typename Frame>
struct ParabolaBasis {
  using frame = Frame;
  using vector_type = VectorType<frame::iN>;
  using arg_vector_type = ArgConstType<vector_type>;

  // using basis_type = fixed_storage_type<bp::iN * 5>;
  using matrix_type = MatrixType<frame::iN, 5>;
  using arg_matrix_type = ArgConstType<matrix_type>;

  // internal variables
  vector_type a, b, c, d, ai, bi, ci;

  // basis:
  //  std::array<double, bp::iN * 5> basis;
  matrix_type basis;

  // spans/views
  ParabolaBasis() = default;
  ~ParabolaBasis() = default;

  constexpr inline void clear() {
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
  constexpr inline void prepare(const arg_vector_type& dx) {
    for (index_type i = frame::i0; i < frame::iN - 1; ++i) {
      a[i] = dx[i] + dx[i + 1];
      b[i] = a[i] + dx[i];
      c[i] = a[i] + dx[i + 1];
      ai[i] = 1. / a[i];
      bi[i] = 1. / b[i];
      ci[i] = 1. / c[i];
    }

    for (index_type i = frame::i0 + 1; i < frame::iN - 2; ++i) {
      //      size_t b_row = i * 5;
      d[i] = 1. / (a[i - 1] + a[i + 1]);
      basis(i, 0) = dx[i] * ai[i] + 2.0 * dx[i + 1] * dx[i] * d[i] * ai[i] *
                                        (a[i - 1] * bi[i] - a[i + 1] * ci[i]);

      basis(i, 1) = -d[i] * dx[i] * a[i - 1] * bi[i];
      basis(i, 2) = d[i] * dx[i + 1] * a[i + 1] * ci[i];

      // check this range - upd: should be fine, one extra term worth not
      // restarting loop;
      d[i] = dx[i] / (a[i - 1] + dx[i + 1]);
      basis(i, 3) = d[i] * b[i - 1] * ai[i];
      basis(i, 4) = d[i] * c[i] * ai[i - 1];
    }
  }
};
template <typename Basis>
struct Parabola {
  using basis_type = Basis;
  using frame = typename basis_type::frame;
  using vector_type = typename basis_type::vector_type;
  using matrix_type = typename basis_type::matrix_type;
  using arg_vector_type = typename basis_type::arg_vector_type;
  using arg_matrix_type = typename basis_type::arg_matrix_type;

  // parabola coeffs
  vector_type da, a6, al, ar;

  // internals
  vector_type diff_an, del_a;
  vector_type diff_amin;
  vector_type tmp1, tmp2, tmp3, onemfl;

  Parabola() = default;
  ~Parabola() = default;

  constexpr inline void clear() {
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

  /*template <typename DerivedV, typename DerivedM>
  constexpr inline void construct(const ArgVectorXpr<DerivedV>& a,
                                  const ArgVectorXpr<DerivedV>& flat,
                                  const ArgMatrixXpr<DerivedM>& pbase) {*/
  /*  template <typename Derived>
    constexpr inline void construct(ArgVectorXpr<Derived>& a,
                                    const vector_type& flat,
                                    const matrix_type& pbase) {*/

  constexpr inline void construct(const arg_vector_type& a,
                                  const arg_vector_type& flat,
                                  const arg_matrix_type& pbase) {
    // clear();

    using Eigen::fix;
    using Eigen::seqN;

    constexpr auto i0 = frame::i0;
    constexpr auto iN = frame::iN;
    constexpr auto nN = iN - i0;

    const auto Nm1 = fix<(iN - i0) - 1>;
    const auto Nm2 = fix<(iN - i0) - 2>;
    const auto Nm3 = fix<(iN - i0) - 3>;
    const auto Nm4 = fix<(iN - i0) - 4>;

    // diff_an(seqN(i0, Nm1)) = a(seqN(i0 + 1, Nm1)) - a(seqN(i0, Nm1));
    // v_range<i0, nN - 1> r1;
    // diff_an(r1) = a(r1.up1()) - a(r1);

    auto r1 = seqN(i0, Nm1);
    diff_an(r1) = a(r1 + 1) - a(r1);

    // v_range<i0 + 1, nN - 2> r2;
    auto r2 = seqN(i0 + 1, Nm2);
    del_a(r2) =
        pbase.col(3)(r2) * diff_an(r2) + pbase.col(4)(r2) * diff_an(r2 - 1);
    /*del_a(seqN(i0 + 1, Nm2)) =
        pbase.col(3)(seqN(i0 + 1, Nm2)) * diff_an(seqN(i0 + 1, Nm2)) +
        pbase.col(4)(seqN(i0 + 1, Nm2)) * diff_an(seqN(i0, Nm2));*/

    diff_amin(r2) = diff_an(r2 - 1).min(diff_an(r2));

    del_a(r2) *= (del_a(r2)).abs().min(two * diff_amin(r2)).sign();
    del_a(r2) = (diff_an(r2 - 1) * diff_an(r2)).select(zero, del_a(r2));
    /*    diff_amin(seqN(i0 + 1, Nm2)) =
            diff_an(seqN(i0, Nm2)).min(diff_an(seqN(i0 + 1, Nm2)));

        del_a(seqN(i0 + 1, Nm2)) *= (del_a(seqN(i0 + 1, Nm2)))
                                        .abs()
                                        .min(two * diff_amin(seqN(i0 + 1, Nm2)))
                                        .sign();

    del_a(seqN(i0 + 1, Nm2)) =
        (diff_an(seqN(i0, Nm2)) * diff_an(seqN(i0 + 1, Nm2)))
            .select(zero, del_a(seqN(i0 + 1, Nm2)));
    */

    //    auto twodam_ = (2.0 * diff_an_.abs());
    //    auto twodamm_ = twodam_.min(twodam_(im1));

    //    del_a = ((diff_an_ * diff_an_(im1)) < zero)
    //                .select(zero, del_a * (del_a.abs().min(twodamm_)).sign());

    // two_diff_an_abs = 2.0 * diff_an.abs();
    // two_diff_an_min = two_diff_an_abs(im1).min(two_diff_an_abs);

    // del_a =
    //    ((diff_an * diff_an(im1)) < zero)
    //        .select(zero, del_a * (del_a.abs().min(two_diff_an_min).sign()));

    // v_range<i0 + 1, nN - 3> r3;
    auto r3 = seqN(i0 + 1, Nm3);
    ar(r3) = a(r3) + pbase.col(0)(r3) * diff_an(r3) +
             pbase.col(1)(r3) * del_a(r3 + 1) + pbase.col(2)(r3) * del_a(r3);

    al(r3 + 1) = ar(r3);
    /*ar(seqN(i0 + 1, Nm3)) =
        a(seqN(i0 + 1, Nm3)) +
        pbase.col(0)(seqN(i0 + 1, Nm3)) * diff_an(seqN(i0 + 1, Nm3)) +
        pbase.col(1)(seqN(i0 + 1, Nm3)) * del_a(seqN(i0 + 2, Nm3)) +
        pbase.col(2)(seqN(i0 + 1, Nm3)) * del_a(seqN(i0 + 1, Nm3));

    al(seqN(i0 + 2, Nm3)) = ar(seqN(i0 + 1, Nm3));*/

    // for (index_type i = frame::i0; i < frame::iN - 1; ++i) {
    // diff_an[i] = a[i + 1] - a[i];
    // two_diff_an_abs[i] = 2.0 * std::abs(diff_an[i]);
    //}

    // for (index_type i = frame::i0 + 1; i < frame::iN - 2; ++i) {
    // two_diff_an_min[i] = std::min(two_diff_an_abs[i - 1],
    // two_diff_an_abs[i]);
    //}

    // for (index_type i = frame::i0 + 1; i < frame::iN - 2; ++i) {
    // del_a[i] = pbase(i, 3) * diff_an[i] + pbase(i, 4) * diff_an[i - 1];

    // del_a[i] = std::copysign(
    // std::min(std::abs(del_a[i]), two_diff_an_min[i]),

    // del_a[i]);

    // if (diff_an[i - 1] * diff_an[i] < 0.0)
    // del_a[i] = 0.0;
    //}

    // for (index_type i = frame::i0 + 1; i < frame::iN - 2; ++i) {
    // ar[i] = a[i] + pbase(i, 0) * diff_an[i] + pbase(i, 1) * del_a[i + 1] +
    // pbase(i, 2) * del_a[i];
    // al[i + 1] = ar[i];
    /*}*/

    // const auto cmi = seqN(i0 + 2, Nm4);
    // v_range<i0 + 2, nN - 4> r4;
    auto r4 = seqN(i0 + 2, Nm4);
    onemfl(r4) = 1.0 - flat(r4);

    ar(r4) = flat(r4) * a(r4) + onemfl(r4) * ar(r4);
    al(r4) = flat(r4) * a(r4) + onemfl(r4) * al(r4);

    //    for (index_type i = frame::i0 + 2; i < frame::iN - 2; ++i) {
    //      double onemfl = 1.0 - flat[i];
    //      ar[i] = flat[i] * a[i] + onemfl * ar[i];
    //      al[i] = flat[i] * a[i] + onemfl * al[i];
    //    }

    da(r4) = ar(r4) - al(r4);
    a6(r4) = (6.0 * (a(r4) - half * (al(r4) + ar(r4))));
    /*    auto tmp1_ = (ar - a) * (a - al);
        auto tmp2_ = da.square();
        auto tmp3_ = da * a6;*/
    tmp1(r4) = (ar(r4) - a(r4)) * (a(r4) - al(r4));
    tmp2(r4) = da(r4).square();
    tmp3(r4) = da(r4) * a6(r4);

    /*    for (index_type i = frame::i0 + 2; i < frame::iN - 2; ++i) {
          //      da[i] = ar[i] - al[i];
          //      a6[i] = 6.0 * (a[i] - half * (al[i] + ar[i]));
          da[i] = t1[i];
          a6[i] = t2[i];

          tmp1[i] = (ar[i] - a[i]) * (a[i] - al[i]);
          tmp2[i] = da[i] * da[i];
          tmp3[i] = da[i] * a6[i];
        }
        */

    for (index_type i = frame::i0 + 2; i < frame::iN - 2; ++i) {
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

    //    for (index_type i = frame::i0 + 2; i < frame::iN - 2; ++i) {
    //      da[i] = ar[i] - al[i];
    //      a6[i] = 6.0 * (a[i] - half * (al[i] + ar[i]));
    //    }
  }

  constexpr inline double rc_left(const index_type i,
                                  const double f1,
                                  const double f2) const {
    return ar[i] - f1 * (da[i] - f2 * a6[i]);
  }

  constexpr inline double rc_right(const index_type i,
                                   const double f1,
                                   const double f2) const {
    return al[i] + f1 * (da[i] + f2 * a6[i]);
  }
};

template <typename ParabolaType, index_type NP>
struct ParabolaSet {
  using parabola_type = ParabolaType;
  using frame = typename parabola_type::frame;
  using vector_type = typename parabola_type::vector_type;
  using arg_vector_type = typename parabola_type::arg_vector_type;
  using arg_matrix_type = typename parabola_type::arg_matrix_type;

  using parabola_basis = typename parabola_type::basis_type;
//  using parabola_array = FluidVars<parabola_type>;
  using parabola_array = std::array<parabola_type, NP>;
  using p_data = std::array<vector_type&, NP>;

  constexpr static index_type NPARA = NP;

  parabola_basis basis;
  parabola_array parabolas;

  ParabolaSet() = default;
  ~ParabolaSet() = default;

  inline void prepare(const arg_vector_type& dx) { basis.prepare(dx); }

  //  template <index_type FN>
  constexpr inline void build(const p_data& fluid_v,
                              const vector_type& flat){
                              //   const std::array<FLUID_VI, FN>& fidx) {) {
    for (index_type i = 0; i < NP; ++i) {
      parabolas[i].construct(fluid_v[i], flat, basis.basis);
    }
  }

  constexpr inline const parabola_type& operator[](const index_type idx) const {
    return parabolas[idx];
  }

  constexpr inline const parabola_array& arr() const { return parabolas; }
};
}  // namespace ppm
