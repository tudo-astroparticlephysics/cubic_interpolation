#include "../Interpolant.h"

#include <boost/math/tools/roots.hpp>

namespace cubic_splines {
namespace detail {
float _find_parameter(std::function<std::tuple<float, float>(float)> func,
                      float x, cubic_splines::Axis const &axis) {
  return boost::math::tools::newton_raphson_iterate(func, x, axis.GetLow(),
                                                    axis.GetHigh(), 10);
}
} // namespace detail
} // namespace cubic_splines
