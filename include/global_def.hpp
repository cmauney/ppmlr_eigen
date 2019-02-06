#pragma once

#include <array>
#include <utility>

#include <Eigen/Dense>

// index type
using index_type = Eigen::Index; // std::ptrdiff_t;

template<index_type N>
using VectorType = typename Eigen::Array<double, N, 1>;

template<index_type N, index_type M>
using MatrixType = typename Eigen::Array<double, N, M>;

template<typename V>
using ArgConstType = Eigen::Ref<const V>;

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
#define EIGEN_NO_AUTOMATIC_RESIZING

// fluid variable IDs
enum
{
  RHO,
  PRS,
  ENT,
  EIN,
  TMP,
  FLT,
  XVL,
  YVL,
  ZVL,
  N_FLUID_COLS
};

// used for convenient Eigen slicing
// note: no bounds checking is done (for now)
// template template parameter to avoid conflicts with (int + int)
template<template<typename...> typename Seq, typename... Ts>
constexpr inline auto
operator+(Seq<Ts...> s, index_type a)
{
  return Eigen::seqN(s.first() + a, s.size());
}

template<template<typename...> typename Seq, typename... Ts>
constexpr inline auto
operator-(Seq<Ts...> s, index_type a)
{
  return Eigen::seqN(s.first() - a, s.size());
}
