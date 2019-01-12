#pragma once

#include <tuple>
#include <utility>
#include <valarray>

#include "global_def.hpp"
#include "ppmlr/ppm_constant.hpp"

namespace ppm
{

constexpr double xwig = 0.0;

constexpr double omega1 = 0.75;
constexpr double omega2 = 5.0;
constexpr double epsilon = 0.33;

constexpr index_type N_STATEV = 3;
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

template <typename Container>
using FluidVars = std::array<Container, N_FLUIDV>;

}  // namespace ppm
