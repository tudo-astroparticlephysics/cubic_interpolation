#pragma once

#include "Axis.h"
#include "InterpolantBuilder.h"

#include <array>
#include <boost/type_traits/make_void.hpp>
#include <memory>
#include <utility>

template <typename T, typename = boost::void_t<>>
struct is_iterable : std::false_type {};

template <typename T>
struct is_iterable<T, boost::void_t<decltype(std::declval<T>().begin(),
                                             std::declval<T>().end())>>
    : std::true_type {};

/**
 * @brief a
 */
template <typename T1, typename T2 = typename T1::Definition>
class Interpolant {
  T2 def;
  T1 inter;

  template <typename T> float back_transform(T trafo, float val) {
    if (trafo)
      val = trafo->back_transform(val);
    return val;
  }

public:
  Interpolant(T2 &&_def, std::string _path, std::string _filename)
      : def(std::move(_def)),
        inter(InterpolantBuilder<T1>().build(def, _path, _filename)){};

  float evaluate(float x) {
    auto x_transformed = def.axis->transform(x);
    auto val = inter.evaluate(x_transformed);
    return back_transform(def.f_trafo.get(), val);
  };

  template <typename T, std::enable_if_t<is_iterable<T>::value, bool> = true>
  float evaluate(T x) {
    auto x_transformed = std::array<float, T1::N>();
    for (size_t i = 0; i < T1::N; ++i)
      x_transformed[i] = def.axis[i]->transform(x[i]);
    auto val = inter.evaluate(x_transformed);
    return back_transform(def.f_trafo.get(), val);
  };
};
