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
 * @brief Utility class for a better handling of interpolants. Take care about
 * transformation of axis building, loading and storage of intepolation tables.
 */
template <typename T1, typename T2 = typename T1::Definition>
class Interpolant {
  T2 def;
  T1 inter;

  template <typename T> float back_transform(T trafo, float val) const {
    if (trafo)
      val = trafo->back_transform(val);
    return val;
  }

public:
  /**
   * @brief To initaialize a Interpolant the corresponding definition to the
   * Interpolant type is required and a path where to store the tables, as well
   * as a filename. If no path and filename is specified, tables will not be
   * stored and a warning thrown because they cann't be reused.
   */
  Interpolant(T2 &&_def, std::string _path = "", std::string _filename = "")
      : def(std::move(_def)),
        inter(InterpolantBuilder<T1>().build(def, _path, _filename)){};

  float evaluate(float x) const {
    auto x_transformed = def.axis->transform(x);
    auto val = inter.evaluate(x_transformed);
    return back_transform(def.f_trafo.get(), val);
  };

  T2 const &GetDefinition() const { return def; };

  template <typename T, std::enable_if_t<is_iterable<T>::value, bool> = true>
  float evaluate(T x) const {
    auto x_transformed = std::array<float, T1::N>();
    for (size_t i = 0; i < T1::N; ++i)
      x_transformed[i] = def.axis[i]->transform(x[i]);
    auto val = inter.evaluate(x_transformed);
    return back_transform(def.f_trafo.get(), val);
  };

  float prime(float x) const {
    auto x_transformed = def.axis->transform(x);
    auto val = inter.prime(x_transformed);
    if (def.f_trafo)
      val *= def.f_trafo->derive(evaluate(x));
    return val / def.axis->derive(x);
  }

  template <typename T, std::enable_if_t<is_iterable<T>::value, bool> = true>
  auto prime(T x) const {
    auto x_transformed = std::array<float, T1::N>();
    for (size_t i = 0; i < T1::N; ++i)
      x_transformed[i] = def.axis[i]->transform(x[i]);
    auto val = inter.prime(x_transformed);
    for (size_t i = 0; i < T1::N; ++i) {
      val[i] *= 1 / def.axis[i]->derive(x[i]);
      if (def.f_trafo)
        val[i] *= def.f_trafo->derive(evaluate(x));
    }
    return val;
  };
};

#include <boost/math/tools/roots.hpp>

template <typename T1, typename T2>
auto FindParameter(T1 const &interpolant, float val, T2 x, size_t n) {
  auto f = [&interpolant, val, &x, n](float xi) {
    x[n] = xi;
    return std::make_tuple(interpolant.evaluate(x) - val,
                           interpolant.prime(x)[n]);
  };
  auto low = interpolant.GetDefinition().axis[n]->GetLow();
  auto high = interpolant.GetDefinition().axis[n]->GetHigh();
  return boost::math::tools::newton_raphson_iterate(f, x[n], low, high, 10);
}

template <typename T1>
auto FindParameter(T1 const &interpolant, float val, float x) {
  auto func = [&interpolant, val](float xi) {
    return std::make_tuple(interpolant.evaluate(xi) - val,
                           interpolant.prime(xi));
  };
  auto low = interpolant.GetDefinition().axis->GetLow();
  auto high = interpolant.GetDefinition().axis->GetHigh();
  return boost::math::tools::newton_raphson_iterate(func, x, low, high, 10);
}
