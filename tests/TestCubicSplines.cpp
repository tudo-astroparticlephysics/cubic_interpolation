
#include "CubicInterpolation/Axis.h"
#include "CubicInterpolation/CubicSplines.h"
#include "CubicInterpolation/Interpolant.h"
#include "gtest/gtest.h"
#include <array>
#include <cmath>
#include <fstream>
#include <random>
#include <vector>

std::random_device rd;
std::mt19937 gen(rd());

/* TEST(CubicSplines, Constructor) { */
/*   auto splineI = cubic_splines::CubicSplines(std::array<double, 5>{1, 2, 3, 4, 5}, 0, 0); */
/*   auto splineII = cubic_splines::CubicSplines(std::vector<double>{1, 2, 3, 4, 5}, 0, 0); */
/*   try { */
/*     auto splineIII = cubic_splines::CubicSplines(std::vector<double>{1, 2}, 0, 0); */
/*   } catch (std::exception const &ex) { */
/*     EXPECT_STREQ( */
/*         ex.what(), */
/*         "Interpolation using a cubic b spline requires at least 3 points.\n"); */
/*   } */
/* } */

TEST(CubicSplines, evaluate_cubic_polynom) {
  size_t N = 5;
  auto low = -10.f;
  auto high = 10.f;
  auto def = cubic_splines::CubicSplines::Definition();
  auto func = [](double x) { return x * x * x - 6 * x * x; };
  def.f = func;
  def.axis = std::make_unique<cubic_splines::LinAxis>(low, high, N);
  auto spline = cubic_splines::Interpolant<cubic_splines::CubicSplines>(std::move(def), "", "");
  std::uniform_real_distribution<double> dis(low, high);
  for (int i = 0; i < 10'000; ++i) {
    auto x = dis(gen);
    EXPECT_NEAR(func(x), spline.evaluate(x),
                std::max(std::abs(func(x)) * 1e-3, 1e-3));
  }
}

TEST(CubicSplines, evaluate_absolute_value) {
  size_t N = 100;
  auto low = -10.f;
  auto high = 10.f;
  auto def = cubic_splines::CubicSplines::Definition();
  auto func = [](double x) { return std::abs(x); };
  def.f = func;
  def.axis = std::make_unique<cubic_splines::LinAxis>(low, high, N);
  auto spline = cubic_splines::Interpolant<cubic_splines::CubicSplines>(std::move(def), "", "");
  std::uniform_real_distribution<double> dis(low, high);
  for (int i = 0; i < 10'000; ++i) {
    auto x = dis(gen);
    EXPECT_NEAR(func(x), spline.evaluate(x),
                std::max(std::abs(func(x)) * 1e-3, 1e-1));
  }
}

TEST(CubicSplines, evaluate_heaviside) {
  size_t N = 100;
  auto low = -0.1f;
  auto high = 0.1f;
  auto def = cubic_splines::CubicSplines::Definition();
  auto func = [](double x) { return x > 0; };
  def.f = func;
  def.axis = std::make_unique<cubic_splines::LinAxis>(low, high, N);
  auto spline = cubic_splines::Interpolant<cubic_splines::CubicSplines>(std::move(def), "", "");
  std::uniform_real_distribution<double> dis(low, high);
  for (int i = 0; i < 10'000; ++i) {
    auto x = dis(gen);
    EXPECT_NEAR(func(x), spline.evaluate(x), 1);
  }
}

TEST(CubicSplines, evaluate_exp_distributed_nodes) {
  size_t N = 100;
  auto low = 1e-3f;
  auto high = 1e1f;
  auto func = [](double x) { return std::exp(x); };
  auto def = cubic_splines::CubicSplines::Definition();
  def.f = func;
  def.axis = std::make_unique<cubic_splines::ExpAxis>(low, high, N);
  auto spline = cubic_splines::Interpolant<cubic_splines::CubicSplines>(std::move(def), "", "");
  std::uniform_real_distribution<double> dis(low, high);

  for (int i = 0; i < 10'000; ++i) {
    auto x = dis(gen);
    EXPECT_NEAR(func(x), spline.evaluate(x), std::abs(func(x)) * 1e-2);
  }
}

TEST(CubicSplines, evaluate_log_func_values_and_axis) {
  size_t N = 3;
  auto low = 1e1f;
  auto high = 1e2f;
  auto def = cubic_splines::CubicSplines::Definition();
  auto func = [](double x) { return std::pow(x, 10); };
  def.f = func;
  def.f_trafo = std::make_unique<cubic_splines::ExpAxis>(1, 0);
  def.axis = std::make_unique<cubic_splines::ExpAxis>(low, high, N);
  auto spline = cubic_splines::Interpolant<cubic_splines::CubicSplines>(std::move(def), "", "");
  std::uniform_real_distribution<double> dis(low, high);

  for (int i = 0; i < 10'000; ++i) {
    auto x = dis(gen);
    EXPECT_NEAR(func(x), spline.evaluate(x), std::abs(func(x)) * 1e-2);
  }
}

TEST(CubicSplines, prime) {
  size_t N = 30;
  auto low = 0.f;
  auto high = 10.f;
  auto def = cubic_splines::CubicSplines::Definition();
  auto func = [](double x) { return x * x + x + 1; };
  auto df_dx = [](double x) { return 2 * x + 1; };
  def.f = func;
  def.axis = std::make_unique<cubic_splines::LinAxis>(low, high, N);
  auto spline = cubic_splines::Interpolant<cubic_splines::CubicSplines>(std::move(def), "", "");
  std::uniform_real_distribution<double> dis(low, high);
  for (int i = 0; i < 10'000; ++i) {
    auto x = dis(gen);
    EXPECT_NEAR(df_dx(x), spline.prime(x), std::abs(func(x)) * 1e-2);
  }
}

TEST(CubicSplines, prime_trafo) {
  size_t N = 30;
  auto low = 1.e-2f;
  auto high = 1.e2f;
  auto def = cubic_splines::CubicSplines::Definition();
  auto func = [](double x) { return x * x + x + 1; };
  auto df_dx = [](double x) { return 2 * x + 1; };
  def.f = func;
  def.f_trafo = std::make_unique<cubic_splines::ExpAxis>(1, 0);
  def.axis = std::make_unique<cubic_splines::ExpAxis>(low, high, N);
  auto spline = cubic_splines::Interpolant<cubic_splines::CubicSplines>(std::move(def), "", "");
  std::uniform_real_distribution<double> dis(low, high);
  for (int i = 0; i < 10'000; ++i) {
    auto x = dis(gen);
    EXPECT_NEAR(df_dx(x), spline.prime(x), std::abs(func(x)) * 1e-2);
    EXPECT_NEAR(df_dx(x), spline.prime(x), std::abs(func(x)) * 1e-2);
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
