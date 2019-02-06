#pragma once

#include <array>
#include <utility>

#include "global_def.hpp"
#include "ppmlr/ppm_constant.hpp"

namespace ppm {

constexpr double xwig = 0.0;

constexpr double omega1 = 0.75;
constexpr double omega2 = 5.0;
constexpr double epsilon = 0.33;

// constexpr index_type N_STATEV = 3;
constexpr index_type NZ_GHOST = 6;

template<index_type N>
constexpr index_type __PPM_N = N + 2 * NZ_GHOST;

template<index_type N>
using __PPMVector = VectorType<__PPM_N<N>>;

template<index_type N, index_type M>
using __PPMMatrix = MatrixType<__PPM_N<N>, M>;

template<index_type N>
using __PPMVecRef = Eigen::Ref<__PPMVector<N>>;

template<index_type N, index_type M>
using __PPMMatRef = Eigen::Ref<__PPMMatrix<N, M>>;

template<index_type N>
using __PPMFluidVars = std::array<__PPMVector<N>, N_FLUID_COLS>;

template<index_type N, index_type NR>
using __PPMFluidRef = std::array<__PPMVecRef<N>, NR>;

enum
{
  QRHO,
  QPRS,
  QXVL,
  NQV
};
enum
{
  URHO,
  UXVL,
  UYVL,
  UZVL,
  UEIN,
  UENT,
  NUV
};

} // namespace ppm
