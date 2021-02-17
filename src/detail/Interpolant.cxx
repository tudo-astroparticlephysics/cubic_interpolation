#include "CubicInterpolation/Interpolant.h"

#include <boost/math/tools/roots.hpp>

namespace cubic_splines {
namespace detail {
double _find_parameter(std::function<double(double)> const &f,
                      std::function<double(double)> const &df, double x,
                      Axis const &axis) {
  auto low = axis.GetLow();
  auto high = axis.GetHigh();
  auto func = [&f, &df](double x) { return std::make_tuple(f(x), df(x)); };
  return boost::math::tools::newton_raphson_iterate(func, x, low, high, 20);
}
} // namespace detail
} // namespace cubic_splines
