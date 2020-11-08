#pragma once

#include "Axis.h"
#include <cmath>
#include <iostream>
#include <memory>
#include <tuple>
#include <type_traits>
#include <vector>

template <typename T> class InterpolationUtility {
public:
  struct Definition {
    std::array<std::unique_ptr<Axis>, T::N> x_trafo;
    std::unique_ptr<Axis> f_trafo = nullptr;
  };

  InterpolationUtility(T &&, Definition &&);

  template <typename T1> float evaluate(T1 x);

private:
  T interpolant;
  Definition def;
};

template <typename T>
InterpolationUtility<T>::InterpolationUtility(T &&_inter, Definition &&_def)
    : interpolant(std::move(_inter)), def(std::move(_def)) {}

struct Interpolant1d {
  static constexpr size_t N = 1;
  float evaluate(size_t, float) { return 1.f; }
};

template <>
template <>
float InterpolationUtility<Interpolant1d>::evaluate<float>(float x) {
  def.x_trafo[0]->transform(x);
  auto node = def.x_trafo[0]->node(x);
  return interpolant.evaluate(node, x);
}

template <typename T>
template <typename T1>
float InterpolationUtility<T>::evaluate(T1 x) {
  auto node = std::array<size_t, T::N>();

  for (size_t i = 0; i < T::N; ++i) {
    def.x_trafo[i]->transform(x[i]);
    node[i] = def.x_trafo[i]->node(x[i]);
  }
  return interpolant.evaluate(node, x);
}
