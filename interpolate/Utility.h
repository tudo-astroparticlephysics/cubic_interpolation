#pragma once

#include "Axis.h"
#include "InterpolantBuilder.h"
#include <cmath>
#include <iostream>
#include <memory>
#include <tuple>
#include <type_traits>
#include <vector>

template <typename T1, typename T2 = typename T1::Definition>
class Interpolant {
  T2 def;
  T1 inter;

public:
  Interpolant(T2 &&_def, std::string _path, std::string _filename)
      : def(std::move(_def)),
        inter(InterpolantBuilder<T1>().build(def, _path, _filename)){};

  template <typename T> float evaluate(T x) {
    auto node = std::array<size_t, T1::N>();
    for (size_t i = 0; i < T1::N; ++i) {
      def.x_trafo[i]->transform(x[i]);
      node[i] = def.x_trafo[i]->node(x[i]);
    }
    return inter.evaluate(node, x);
  }
};

/* template <> */
/* template <> */
/* float InterpolationUtility<Interpolant1d>::evaluate<float>(float x) { */
/*   def.x_trafo[0]->transform(x); */
/*   auto node = def.x_trafo[0]->node(x); */
/*   return interpolant.evaluate(node, x); */
/* } */

/* template <typename T> */
/* template <typename T1> */
/* float InterpolationUtility<T>::evaluate(T1 x) { */
/*   auto node = std::array<size_t, T::N>(); */

/*   for (size_t i = 0; i < T::N; ++i) { */
/*     def.x_trafo[i]->transform(x[i]); */
/*     node[i] = def.x_trafo[i]->node(x[i]); */
/*   } */
/*   return interpolant.evaluate(node, x); */
/* } */
