#include <cstdlib>
#include <iostream>
#include <memory>

#include "CubicInterpolation/version.h"
#include "CubicInterpolation/Axis.h"
#include "CubicInterpolation/CubicSplines.h"
#include "CubicInterpolation/Interpolant.h"

double func(double x) { return x * x + x; }

int main(int argc, char *argv[]) {

  if (argc != 2) {
    std::cout << "Cubic Interpolation (version " << getCubicInterpolationVersion()
              << ") \nUsage: 1dim_example POINT\n\n  POINT     point to evaluate "
              << "the function\n\n";
    return 0;
  }

  auto def = cubic_splines::CubicSplines<double>::Definition();
  def.f = func;
  def.axis = std::make_unique<cubic_splines::LinAxis<double>>(-2.f, 2.f, (size_t)10);

  auto inter = cubic_splines::Interpolant<cubic_splines::CubicSplines<double>>(
      std::move(def), "", "");

  auto point = std::atof(argv[1]);
  auto res = inter.evaluate(point);

  std::cout << "func(" << point << "): " << func(point) << std::endl;
  std::cout << "f(" << point << "): " << res << std::endl;

  return 0;
}
