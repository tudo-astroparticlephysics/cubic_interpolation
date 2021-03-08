#include "CubicInterpolation/Axis.h"
#include "CubicInterpolation/BicubicSplines.h"
#include "CubicInterpolation/CubicSplines.h"
#include "CubicInterpolation/FindParameter.hpp"
#include "CubicInterpolation/Interpolant.h"
#include "gtest/gtest.h"
#include <array>
#include <iostream>
#include <random>

std::random_device rd;
std::mt19937 gen(rd());

TEST(find_parameter, CubicSplines) {
  constexpr static size_t N = 100;
  auto func = [](double x) { return x * x * x; };
  auto def = cubic_splines::CubicSplines<double>::Definition();
  def.f = func;
  auto xaxis = cubic_splines::LinAxis<double>(-1, 1, N);
  def.axis = std::make_unique<cubic_splines::LinAxis<double>>(xaxis);
  auto spline = cubic_splines::Interpolant<cubic_splines::CubicSplines<double>>(
      std::move(def), "", "");
  std::uniform_real_distribution<double> dis(0, N - 1);
  for (int i = 0; i < 10'000; ++i) {
    auto x = xaxis.back_transform(dis(gen));
    auto f = func(x);
    // search y falue for given f, x value.
    auto guess = cubic_splines::ParameterGuess<double>{.x = 0.5f};
    auto x_guess = cubic_splines::find_parameter(spline, f, guess);
    EXPECT_NEAR(x, x_guess, std::max(std::abs(x) * 1e-2, 1e-3));
  }
}

TEST(find_parameter, BicubicSplines) {
  constexpr static size_t N = 100;
  auto func = [](double x1, double x2) { return x1 * x1 + x2; };
  auto def = cubic_splines::BicubicSplines<double>::Definition();
  def.f = func;
  def.f_trafo = std::make_unique<cubic_splines::ExpAxis<double>>(1, 0);
  auto xaxis = cubic_splines::LinAxis<double>(-1, 1, N);
  auto yaxis = cubic_splines::ExpAxis<double>(1.e0, 1e2, N);
  def.axis[0] = std::make_unique<cubic_splines::LinAxis<double>>(xaxis);
  def.axis[1] = std::make_unique<cubic_splines::ExpAxis<double>>(yaxis);
  auto spline = cubic_splines::Interpolant<cubic_splines::BicubicSplines<double>>(
      std::move(def), "", "");
  std::uniform_real_distribution<double> dis(0, N - 1);
  for (int i = 0; i < 10'000; ++i) {
    auto x = xaxis.back_transform(dis(gen));
    auto y = yaxis.back_transform(dis(gen));
    auto f = func(x, y);
    // search y value for given f, x value.
    auto guess = cubic_splines::ParameterGuess<std::array<double, 2>>{
        .x = std::array<double, 2>{x, 0.5f}, .n = 1, .lower = NAN, .upper = NAN};
    auto y_guess = cubic_splines::find_parameter(spline, f, guess);
    EXPECT_NEAR(y, y_guess, std::max(std::abs(y) * 1e-2, 1e-3));
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
