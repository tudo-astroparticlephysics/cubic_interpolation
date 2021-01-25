#include <iostream>
#include <memory>
#include <string>

#include "CubicInterpolation/Axis.h"
#include "CubicInterpolation/BicubicSplines.h"
#include "CubicInterpolation/Interpolant.h"

double func(double x_1, double x_2) { return x_1 * x_1 + x_2 * x_2 + x_1 * x_2; }

int main(int argc, char *argv[]) {
  if (argc != 3) {
    std::cout << "Usage: 1dim_example POINT_1 POINT_2\n"
              << "\n"
              << "  POINT_1     first variable to evaluate the function\n"
              << "  POINT_2     second variable to evaluate the function\n"
              << std::endl;
    return 0;
  }

  auto def = cubic_splines::BicubicSplines<double>::Definition();
  def.f = func;
  def.axis[0] = std::make_unique<cubic_splines::LinAxis<double>>(-1., 1., (size_t)11);
  def.axis[1] = std::make_unique<cubic_splines::LinAxis<double>>(-1., 1., (size_t)11);

  auto inter =
      cubic_splines::Interpolant<cubic_splines::BicubicSplines<double>>(std::move(def), "/tmp", "interpolant.txt");

  auto point = std::array<double, 2>{std::stof(argv[1]), std::stof(argv[2])};
  auto res = inter.evaluate(point);

  std::cout << "f(" << point[0] << ", " << point[1] << "): " << res
            << std::endl;
  return 0;
}
