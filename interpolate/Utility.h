#pragma once

#include "Axis.h"
#include "InterpolantBuilder.h"
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <type_traits>

template <typename T, typename = std::void_t<>>
struct is_iterable : std::false_type {};

template <typename T>
struct is_iterable<T, std::void_t<decltype(std::declval<T>().begin(),
                                           std::declval<T>().end())>>
    : std::true_type {};

template <typename T1, typename T2 = typename T1::Definition>
class Interpolant {
  T2 def;
  T1 inter;

public:
  Interpolant(T2 &&_def, std::string _path, std::string _filename)
      : def(std::move(_def)),
        inter(InterpolantBuilder<T1>().build(def, _path, _filename)){};

  template <typename T,
            std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
  float evaluate(T x) {
    auto x_rel = def.axis->transform(x);
    auto node = def.axis->node(x_rel);
    return inter.evaluate(node, x_rel);
  };

  template <typename T, std::enable_if_t<is_iterable<T>::value, bool> = true>
  float evaluate(T x) {
    auto node = std::array<size_t, T1::N>();
    auto x_rel = std::array<float, T1::N>();
    for (size_t i = 0; i < T1::N; ++i) {
      x_rel[i] = def.axis[i]->transform(x[i]);
      node[i] = def.axis[i]->node(x_rel[i]);
    }
    return inter.evaluate(node, x_rel);
  };
};
