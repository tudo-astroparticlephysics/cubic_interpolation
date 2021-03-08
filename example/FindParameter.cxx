#include <cmath>
#include <cstdlib>
#include <array>
#include <iostream>
#include <memory>

#include "CubicInterpolation/Axis.h"
#include "CubicInterpolation/BicubicSplines.h"
#include "CubicInterpolation/CubicSplines.h"
#include "CubicInterpolation/FindParameter.hpp"
#include "CubicInterpolation/Interpolant.h"

double func1dim(double x) { return x * x + x; }
double func2dim(double x, double y) { return x * x + y * y; }

int main() {
  auto def1d = cubic_splines::CubicSplines<double>::Definition();
  def1d.f = func1dim;
  def1d.axis = std::make_unique<cubic_splines::LinAxis<double>>(-2.f, 2.f, (size_t)10);
  auto inter1d = cubic_splines::Interpolant<cubic_splines::CubicSplines<double>>(
      std::move(def1d), "", "");
  auto function_value = 1.f;
  auto guess = cubic_splines::ParameterGuess<double>{.x = 1.f};

  auto x = find_parameter(inter1d, function_value, guess);
  std::cout << "f(x) = " << function_value << " -> x = " << x << std::endl;

  auto def2d = cubic_splines::BicubicSplines<double>::Definition();
  def2d.f = func2dim;
  def2d.axis[0] =
      std::make_unique<cubic_splines::LinAxis<double>>(1.e-5f, 1.e1f, (size_t)10);
  def2d.axis[1] =
      std::make_unique<cubic_splines::LinAxis<double>>(1.e-5f, 1.e1f, (size_t)10);
  auto inter2d = cubic_splines::Interpolant<cubic_splines::BicubicSplines<double>>(
      std::move(def2d), "", "");

  function_value = 42.f;
  auto guess_2d = cubic_splines::ParameterGuess<std::array<double, 2>>{
      .x = std::array<double, 2>{0.5f, 3.f}, .n = 1};
  auto y = cubic_splines::find_parameter(inter2d, function_value, guess_2d);
  std::cout << "f(x=" << guess_2d.x[0] << ",y) = " << function_value << " -> y = " << y
            << std::endl;
  return 0;
}
