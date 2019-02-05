#pragma once

#include <tuple>
#include <utility>
#include <array>

#include "global_def.hpp"
#include "ppmlr/ppm_constant.hpp"

namespace ppm
{

constexpr double xwig = 0.0;

constexpr double omega1 = 0.75;
constexpr double omega2 = 5.0;
constexpr double epsilon = 0.33;

//constexpr index_type N_STATEV = 3;
constexpr index_type NZ_GHOST = 6;


// constexpr size_t GRID_NX = 100;
// constexpr size_t GRID_NY = 100;

template <index_type N>
struct VectorFrame {
  constexpr static index_type i0 = 0;
  constexpr static index_type iN = N + 2 * NZ_GHOST;
  constexpr static index_type j0 = NZ_GHOST;
  constexpr static index_type jM = N + NZ_GHOST;
};

template <index_type N>
using PPMGridSet = VectorType<VectorFrame<N>::iN>;

template <index_type N>
using PPMVarSet = MatrixType<VectorFrame<N>::iN, N_FLUID_COLS>;

enum {QRHO, QPRS, QXVL, NQV};
enum {URHO, UXVL, UYVL, UZVL, UEIN, UENT, NUV};

template<index_type N, index_type NV>
using PPMSub = std::array <Eigen::Ref
                            < 
                              VectorType<VectorFrame<N>::iN> 
                            >, 
                          NV
                          >;


//template <typename Container>
//using FluidVars = std::array<Container, N_FLUIDV>;

}  // namespace ppm
