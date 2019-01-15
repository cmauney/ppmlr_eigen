#pragma once

#include <array>
#include <utility>

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

// index type
using index_type = Eigen::Index;  // std::ptrdiff_t;

template <index_type N>
using VectorType = typename Eigen::Array<double, N, 1>;

template <index_type N, index_type M>
using MatrixType = typename Eigen::Array<double, N, M>;

template <typename V>
using ArgConstType = Eigen::Ref<const V>;

template <typename Derived>
using ArgVectorXpr = Eigen::ArrayBase<Derived>;

template <typename Derived>
using ArgMatrixXpr = Eigen::MatrixBase<Derived>;

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
#define EIGEN_NO_AUTOMATIC_RESIZING

// fluid variable IDs
enum { RHO, ENT, EIN, FLT, N_FLUIDS };
enum { VEL, N_FLUIDV };


// IDs of the primative variables (technically prs is not primative but easier
// to use. these define the states of the riemann problem during ppm
//constexpr std::array<FLUID_VI, 3> PRIMATIVE_VARS{RHO, PRS, XVL};

// conserved quantities. note these are 'specific', i.e. per-unit-mass.
//constexpr std::array<FLUID_VI, 6> CONSERVATIVE_VARS{RHO, XVL, YVL,
//                                                    ZVL, EIN, ENT};

// total number of advected quantities
//constexpr index_type N_FLUIDV = 8;

// used for convenient Eigen slicing
// note: no bounds checking is done (for now)
// template template parameter to avoid conflicts with (int + int)
template <template <typename...> typename Seq, typename... Ts>
constexpr inline auto operator+(Seq<Ts...> s, index_type a) {
  return Eigen::seqN(s.first() + a, s.size());
}

template <template <typename...> typename Seq, typename... Ts>
constexpr inline auto operator-(Seq<Ts...> s, index_type a) {
  return Eigen::seqN(s.first() - a, s.size());
}

template<index_type ...NXs>
struct GridSet
{
  constexpr static index_type NDIM = sizeof...(NXs);
  std::array<VectorType<NXs>..., NDIM> xe, xc, dx;
};

template<index_type ...NXs>
struct VarSet
{
  using scalar_field = Eigen::TensorFixedSize<double, Eigen::Sizes<NXs...>>;
  
  using vector_field = Eigen::TensorFixedSize<double, Eigen::Sizes<NXs...>>;

  std::array<scalar_field, N_FLUIDS> fluid_sca;
  std::array<vector_field, N_FLUIDV> fluid_vec;
  
};
