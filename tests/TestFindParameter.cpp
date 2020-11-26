#include "Axis.h"
#include "BicubicSplines.h"
#include "Utility.h"
#include "gtest/gtest.h"
#include <iostream>
#include <array>
#include <random>

std::random_device rd;
std::mt19937 gen(rd());

TEST(FindParameter, CubicSplines) {
  constexpr static size_t N = 100;
  auto func = [](float x) { return x * x * x ; };
  auto def = CubicSplines::Definition();
  def.f = func;
  auto xaxis = LinAxis(-1, 1, N);
  def.axis = std::make_unique<LinAxis>(xaxis);
  auto spline = Interpolant<CubicSplines>(std::move(def), "", "");
  std::uniform_real_distribution<float> dis(0, N-1);
  for (int i = 0; i < 10'000; ++i) {
    auto x = xaxis.back_transform(dis(gen));
    auto f = func(x);
    // search y falue for given f, x value.
    auto x_guess = FindParameter(spline, f, 0.5f);
    EXPECT_NEAR(x, x_guess, std::max(std::abs(x) * 1e-2, 1e-3));
  }
}

TEST(FindParameter, BicubicSplines) {
  constexpr static size_t N = 100;
  auto func = [](float x1, float x2) { return x1 * x1 + x2; };
  auto def = BicubicSplines::Definition();
  def.f = func;
  def.f_trafo = std::make_unique<ExpAxis>(1, 0);
  auto xaxis = LinAxis(-1, 1, N);
  auto yaxis = ExpAxis(1.e0, 1e2, N);
  def.axis[0] = std::make_unique<LinAxis>(xaxis);
  def.axis[1] = std::make_unique<ExpAxis>(yaxis);
  auto spline = Interpolant<BicubicSplines>(std::move(def), "", "");
  std::uniform_real_distribution<float> dis(0, N-1);
  for (int i = 0; i < 10'000; ++i) {
    auto x = xaxis.back_transform(dis(gen));
    auto y = yaxis.back_transform(dis(gen));
    auto f = func(x, y);
    // search y falue for given f, x value.
    auto y_guess = FindParameter(spline, f, std::array<float, 2>{x, 0.5f}, 1);
    EXPECT_NEAR(y, y_guess, std::max(std::abs(y) * 1e-2, 1e-3));
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
