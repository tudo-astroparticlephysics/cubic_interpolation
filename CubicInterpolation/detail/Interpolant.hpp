#include "../Axis.h"

#include <functional>
#include <iostream>

namespace cubic_splines {
namespace detail {

typedef Axis<double> Axis;

template <typename T, std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
auto transform(Axis const &axis, T x) {
  return axis.transform(x);
}

template <typename T1, typename T2> auto transform(T1 const &axis, T2 x) {
  for (auto i = 0u; i < axis.size(); ++i)
    x[i] = axis[i]->transform(x[i]);
  return x;
}

template <typename T> auto back_transform(T trafo, double val) {
  if (trafo)
    val = trafo->back_transform(val);
  return val;
}

template <typename T>
auto back_transform_prime(T trafo, Axis const &axis, double f, double df, double x) {
  if (trafo)
    df *= f;
  df *= axis.derive(x);
  return df;
}

template <typename T1, typename T2, typename T3>
auto back_transform_prime(T1 trafo, T2 const &axis, double f, T3 df, T3 x) {
  for (auto i = 0u; i < axis.size(); ++i)
    df[i] = back_transform_prime(trafo, *axis[i], f, df[i], x[i]);
  return df;
}

double _find_parameter(std::function<double(double)> const &f,
                       std::function<double(double)> const &df, double x,
                       Axis const &axis);

template <typename T1> auto find_parameter(T1 const &inter, double val, double x) {
  auto f = [&inter, val](double xi) { return inter.evaluate(xi) - val; };
  auto df = [&inter](double xi) { return inter.prime(xi); };
  return _find_parameter(f, df, x, *inter.GetDefinition().axis);
};

template <typename T> auto updated_val(T &x, size_t n, double xi) {
  x[n] = xi;
  return x;
}

template <typename T1, typename T2>
auto find_parameter(T1 const &inter, double val, T2 x, size_t n) {
  auto f = [&inter, val, &x, n](double xi) {
    return inter.evaluate(updated_val(x, n, xi)) - val;
  };
  auto df = [&inter, &x, n](double xi) { return inter.prime(updated_val(x, n, xi))[n]; };
  return _find_parameter(f, df, x[n], *inter.GetDefinition().axis[n]);
}
} // namespace detail
} // namespace cubic_splines
