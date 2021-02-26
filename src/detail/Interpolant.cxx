#include "CubicInterpolation/Interpolant.h"

#include <boost/math/tools/roots.hpp>
#include <cassert>

namespace cubic_splines {
namespace detail {
double _find_parameter(std::function<double(double)> const &f,
                       std::function<double(double)> const &df, Axis const &axis,
                       double initial_guess) {
  auto low = axis.GetLow();
  auto high = axis.GetHigh();
  assert((low == high) && "find_parameter call is not well defined");
  auto func = [&f, &df](double x) { return std::make_tuple(f(x), df(x)); };

  if (std::isnan(initial_guess)) {
    auto tol = [](double min, double max) {
      return (std::abs(max - min) < 1.e-2) ? true : false;
    };
    std::tie(low, high) = boost::math::tools::bisect(f, low, high, tol);
    initial_guess = (low + high) / 2;
  }

  return boost::math::tools::newton_raphson_iterate(func, initial_guess, low, high, 20);
}
} // namespace detail
} // namespace cubic_splines
