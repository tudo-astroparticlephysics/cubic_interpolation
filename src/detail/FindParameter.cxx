
#include "CubicInterpolation/FindParameter.hpp"

#include <boost/math/tools/roots.hpp>
#include <cassert>

namespace cubic_splines {
namespace detail {
double _find_parameter(std::function<double(double)> const &f,
                        std::function<double(double)> const &df, double x_guess,
                        double lower, double upper) {
  assert(((void)"find_parameter call is not well defined", lower < upper));
  auto func = [&f, &df](double x) { return std::make_tuple(f(x), df(x)); };

  if (std::isnan(x_guess)) {
    auto bisec_tolerance = (upper - lower) * 1e-2;
    auto tol = [&bisec_tolerance](double min, double max) {
      return (std::abs(max - min) < bisec_tolerance) ? true : false;
    };
    std::tie(lower, upper) = boost::math::tools::bisect(f, lower, upper, tol);
    x_guess = (lower + upper) / 2;
    if (lower == upper)
      return lower;
  }

  return boost::math::tools::newton_raphson_iterate(func, x_guess, lower, upper, 20);
}
} // namespace detail
} // namespace cubic_splines
