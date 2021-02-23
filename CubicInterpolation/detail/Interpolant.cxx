#include "../Interpolant.h"

#include <boost/math/tools/roots.hpp>

namespace cubic_splines {
namespace detail {
double _find_parameter(std::function<double(double)> const &f,
                      std::function<double(double)> const &df, double x,
                      Axis const &axis) {
  auto low = axis.GetLow();
  auto high = axis.GetHigh();
  auto func = [&f, &df](double x) { return std::make_tuple(f(x), df(x)); };

  auto tol = [](double min, double max) { return (std::abs(max - min) < 1.e-2) ? true : false; };
  std::tie(low, high) = boost::math::tools::bisect(f, low, high, tol); 

  return boost::math::tools::newton_raphson_iterate(func, (low + high)/2, low, high, 20);
}
} // namespace detail
} // namespace cubic_splines
