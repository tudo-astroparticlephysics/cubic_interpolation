#include "Axis.h"
#include "BicubicSplines.h"
#include "Utility.h"
#include "gtest/gtest.h"
#include <array>
#include <cmath>
#include <fstream>
#include <random>
#include <vector>

std::random_device rd;
std::mt19937 gen(rd());

TEST(BicubicSplines, Constructor) {
  auto def = BicubicSplines::Definition();
  def.f = [](float x1, float x2) { return x1 * x1 + x2 * x2; };
  def.axis[0] = std::make_unique<LinAxis>(-1, 1, (size_t)10);
  def.axis[1] = std::make_unique<LinAxis>(-1, 1, (size_t)10);

  auto spline = Interpolant<BicubicSplines>(std::move(def), "", "");
};

TEST(BicubicSplines, evaluate_cubic_polynom) {
  size_t N = 30;
  auto low = -10.f;
  auto high = 10.f;
  auto def = BicubicSplines::Definition();
  auto func = [](float x1, float x2) { return x1 * x1 + x2 * x2; };
  def.f = func;
  def.axis[0] = std::make_unique<LinAxis>(low, high, N);
  def.axis[1] = std::make_unique<LinAxis>(low, high, N);
  auto spline = Interpolant<BicubicSplines>(std::move(def), "", "");
  std::uniform_real_distribution<float> dis(low, high);
  for (int i = 0; i < 10'000; ++i) {
    auto x = dis(gen);
    auto y = dis(gen);
    EXPECT_NEAR(func(x, y), spline.evaluate(std::array<float, 2>{x, y}),
                std::max(std::abs(func(x, y)) * 1e-3, 1e-3));
  }
};

TEST(BicubicSplines, evaluate_absolute_value) {
  size_t N = 30;
  auto low = -10.f;
  auto high = 10.f;
  auto def = BicubicSplines::Definition();
  auto func = [](float x1, float x2) { return std::abs(x1) + std::abs(x2); };
  def.f = func;
  def.axis[0] = std::make_unique<LinAxis>(low, high, N);
  def.axis[1] = std::make_unique<LinAxis>(low, high, N);
  auto spline = Interpolant<BicubicSplines>(std::move(def), "", "");
  std::uniform_real_distribution<float> dis(low, high);
  for (int i = 0; i < 10'000; ++i) {
    auto x = dis(gen);
    auto y = dis(gen);
    EXPECT_NEAR(func(x, y), spline.evaluate(std::array<float, 2>{x, y}),
                std::max(std::abs(func(x, y)) * 1e-3, 5e-1));
  }
}

TEST(BicubicSplines, evaluate_heaviside) {
  size_t N = 30;
  auto low = -0.1f;
  auto high = 0.1f;
  auto def = BicubicSplines::Definition();
  auto func = [](float x1, float x2) { return ((x1 > 0) && (x2 > 0)); };
  def.f = func;
  def.axis[0] = std::make_unique<LinAxis>(low, high, N);
  def.axis[1] = std::make_unique<LinAxis>(low, high, N);
  auto spline = Interpolant<BicubicSplines>(std::move(def), "", "");
  std::uniform_real_distribution<float> dis(low, high);
  for (int i = 0; i < 10'000; ++i) {
    auto x = dis(gen);
    auto y = dis(gen);
    EXPECT_NEAR(func(x, y), spline.evaluate(std::array<float, 2>{x, y}), 1);
  }
}

TEST(BicubicSplines, evaluate_exp_distributed_nodes) {
  size_t N = 30;
  auto low = 1e-5f;
  auto high = 1e0f;
  auto def = BicubicSplines::Definition();
  auto func = [](float x1, float x2) { return std::exp(x1 * x1 + x2 * x2); };
  def.f = func;
  def.axis[0] = std::make_unique<ExpAxis>(low, high, N);
  def.axis[1] = std::make_unique<ExpAxis>(low, high, N);
  auto spline = Interpolant<BicubicSplines>(std::move(def), "", "");
  std::uniform_real_distribution<float> dis(low, high);
  for (int i = 0; i < 10'000; ++i) {
    auto x = dis(gen);
    auto y = dis(gen);
    EXPECT_NEAR(func(x, y), spline.evaluate(std::array<float, 2>{x, y}), 1);
  }
}

TEST(BicubicSplines, evaluate_log_func_values_and_axis) {
  size_t N = 10;
  auto low = 1e0f;
  auto high = 1e1f;
  auto def = BicubicSplines::Definition();
  auto func = [](float x1, float x2) {
    return std::pow(x1, 10) * std::pow(x2, 10);
  };
  def.f = func;
  def.f_trafo = std::make_unique<ExpAxis>(1, 0);
  def.axis[0] = std::make_unique<ExpAxis>(low, high, N);
  def.axis[1] = std::make_unique<ExpAxis>(low, high, N);
  auto spline = Interpolant<BicubicSplines>(std::move(def), "", "");
  std::uniform_real_distribution<float> dis(low, high);
  for (int i = 0; i < 10; ++i) {
    auto x = dis(gen);
    auto y = dis(gen);
    EXPECT_NEAR(func(x, y), spline.evaluate(std::array<float, 2>{x, y}),
                std::abs(func(x, y)) * 1e-2);
  }
}

TEST(BicubicSplines, prime) {
  size_t N = 11;
  auto low = 1.e0f;
  auto high = 1.e1f;
  auto def = BicubicSplines::Definition();
  auto func = [](float x1, float x2) { return x1 * x1 + x2 * x2 + x1 * x2; };
  auto df_dx = [](float x1, float x2) {
    return std::array<float, 2>{2 * x1 + x2, 2 * x2 + x1};
  };
  def.f = func;
  def.axis[0] = std::make_unique<LinAxis>(low, high, N);
  def.axis[1] = std::make_unique<LinAxis>(low, high, N);
  auto splineI = Interpolant<BicubicSplines>(std::move(def), "", "");
  std::uniform_real_distribution<float> dis(low, high);
  for (int i = 0; i < 10'000; ++i) {
    auto x = dis(gen);
    auto y = dis(gen);
    EXPECT_NEAR(func(x, y), splineI.evaluate(std::array<float, 2>{x, y}),
                std::abs(func(x, y)) * 1e-2);
    auto f = df_dx(x, y);
    auto sp = splineI.prime(std::array<float, 2>{x, y});
    EXPECT_NEAR(f[0], sp[0], std::abs(f[0]) * 1e-2);
    EXPECT_NEAR(f[1], sp[1], std::abs(f[1]) * 1e-2);
  }

  def = BicubicSplines::Definition();
  def.f = func;
  def.axis[0] = std::make_unique<ExpAxis>(low, high, N);
  def.axis[1] = std::make_unique<ExpAxis>(low, high, N);
  auto splineII = Interpolant<BicubicSplines>(std::move(def), "", "");
  for (int i = 0; i < 10'000; ++i) {
    auto x = dis(gen);
    auto y = dis(gen);
    EXPECT_NEAR(func(x, y), splineII.evaluate(std::array<float, 2>{x, y}),
                std::abs(func(x, y)) * 1e-2);
    auto f = df_dx(x, y);
    auto sp = splineII.prime(std::array<float, 2>{x, y});
    EXPECT_NEAR(f[0], sp[0], std::abs(f[0]) * 1e-2);
    EXPECT_NEAR(f[1], sp[1], std::abs(f[1]) * 1e-2);
  }
}

TEST(BicubicSplines, prime_trafo) {
  size_t N = 11;
  auto low = 1.e0f;
  auto high = 1.e1f;
  auto def = BicubicSplines::Definition();
  auto func = [](float x1, float x2) { return x1 * x1 + x2 * x2 + x1 * x2; };
  auto df_dx = [](float x1, float x2) {
    return std::array<float, 2>{2 * x1 + x2, 2 * x2 + x1};
  };
  def.f = func;
  def.f_trafo = std::make_unique<ExpAxis>(1, 0);
  def.axis[0] = std::make_unique<ExpAxis>(low, high, N);
  def.axis[1] = std::make_unique<LinAxis>(low, high, N);
  auto spline = Interpolant<BicubicSplines>(std::move(def), "", "");
  std::uniform_real_distribution<float> dis(low, high);
  for (int i = 0; i < 10'000; ++i) {
    auto x = dis(gen);
    auto y = dis(gen);
    EXPECT_NEAR(func(x, y), spline.evaluate(std::array<float, 2>{x, y}),
                std::abs(func(x, y)) * 1e-2);
    auto f = df_dx(x, y);
    auto sp = spline.prime(std::array<float, 2>{x, y});
    EXPECT_NEAR(f[0], sp[0], std::abs(f[0]) * 1e-2);
    EXPECT_NEAR(f[1], sp[1], std::abs(f[1]) * 1e-2);
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
