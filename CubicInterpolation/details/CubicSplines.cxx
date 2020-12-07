#include "../CubicSplines.h"

namespace cubic_splines {

float CubicSplines::evaluate(float x) const { return spline(x); };

float CubicSplines::prime(float x) const { return spline.prime(x); };

float CubicSplines::double_prime(float x) const {
  return spline.double_prime(x);
};

CubicSplines CubicSplines::Data::build() const {
  return CubicSplines(y, lower_lim_derivate, upper_lim_derivate);
};
} // namespace cubic_splines
