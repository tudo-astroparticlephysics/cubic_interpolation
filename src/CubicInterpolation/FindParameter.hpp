#pragma once

#include "CubicInterpolation/Axis.h"

#include <cmath>
#include <functional>
#include <type_traits>

namespace cubic_splines {
template <typename T> struct ParameterGuess {
  T x;
  size_t n = 0;
  double lower = NAN;
  double upper = NAN;
};

namespace detail {

template <typename T, typename = void> struct is_iterable : std::false_type {};

template <typename T>
struct is_iterable<T, std::void_t<decltype(std::declval<T>().begin()),
                                  decltype(std::declval<T>().end())>> : std::true_type {};

template <typename T1, typename T2, std::enable_if_t<is_iterable<T1>::value, bool> = true>
auto PopulateParameterGuess(cubic_splines::ParameterGuess<T1> &guess, T2 const &ax) {
  if (std::isnan(guess.lower))
    guess.lower = ax[guess.n]->GetLow();
  if (std::isnan(guess.upper))
    guess.upper = ax[guess.n]->GetHigh();
}

template <typename T1, typename T2,
          std::enable_if_t<std::is_floating_point<T1>::value, bool> = true>
auto PopulateParameterGuess(cubic_splines::ParameterGuess<T1> &guess, T2 const &ax) {
  if (std::isnan(guess.lower))
    guess.lower = ax.GetLow();
  if (std::isnan(guess.upper))
    guess.upper = ax.GetHigh();
}

template <typename T, std::enable_if_t<is_iterable<T>::value, bool> = true>
auto updated_val(T &x, size_t n, double xi) {
  x[n] = xi;
  return x;
}

template <typename T, std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
auto updated_val(T &x, size_t n, double xi) {
  x = xi;
  return x;
}

double _find_parameter(std::function<double(double)> const &f,
                       std::function<double(double)> const &df, double x_guess,
                       double lower, double upper);

template <typename T, typename T2, std::enable_if_t<is_iterable<T>::value, bool> = true>
auto find_parameter(std::function<double(double)> f, T2 df,
                    cubic_splines::ParameterGuess<T> guess) {
  return _find_parameter(
      f, [df, n = guess.n](double x) { return df(x)[n]; }, guess.x[guess.n], guess.lower,
      guess.upper);
};

template <typename T, std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
auto find_parameter(std::function<double(double)> f, std::function<double(double)> df,
                    cubic_splines::ParameterGuess<T> guess) {
  return _find_parameter(f, df, guess.x, guess.lower, guess.upper);
};
} // namespace detail

/**
 * @brief Inverts the interpolation evaluation and searches the missing axis-
 * for the given function-value. If a n > 1 dimensional interpolant is passed,
 * the parameterspace must be narrowed. The function values except the searched
 * dimension must be specified and a guess must be given in place of the
 * container.
 *
 * @param interpolant interpolant object required for calculation
 * @param val function value
 * @param x axis values including guess of the parameter to estimate. If
 * parameter to estimate is set to `nan`, the parameterspace for the missing
 * axis value will be shrinked automatically the searched xaxis value.
 * @param low lower limit of parameter to estimate. If set to `nan`, lowest
 * possible boundary (given by axis) will be assumed
 * @param high upper limit of parameter to estimate. If set to `nan`,
 * highest possible boundary (given by axis) will be used
 * @param n (only if dim > 1) n-th axis to estimate the value
 */
template <typename T1, typename T2>
auto find_parameter(T1 const &inter, double val, ParameterGuess<T2> guess) {
  detail::PopulateParameterGuess(guess, inter.GetDefinition().GetAxis());
  auto f = [&inter, val, &guess](double xi) {
    return inter.evaluate(detail::updated_val(guess.x, guess.n, xi)) - val;
  };
  auto df = [&inter, &guess](double xi) {
    return inter.prime(detail::updated_val(guess.x, guess.n, xi));
  };
  return detail::find_parameter(f, df, guess);
}
} // namespace cubic_splines
