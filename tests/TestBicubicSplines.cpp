#include "CubicInterpolation/Axis.h"
#include "CubicInterpolation/BicubicSplines.h"
#include "CubicInterpolation/Interpolant.h"
#include "gtest/gtest.h"
#include <array>
#include <cmath>
#include <fstream>
#include <random>
#include <vector>

std::random_device rd;
std::mt19937 gen(rd());

using spline_t = cubic_splines::BicubicSplines<double>;
using spline_def_t = cubic_splines::BicubicSplines<double>::Definition;

class BicubicSplinesTest : public spline_t {
public:
  BicubicSplinesTest(spline_def_t const &_def) : spline_t(_def) {}

  template <typename T1>
  inline auto _prime(Definition const &def, T1 func, unsigned int n1, unsigned int n2) {
    return spline_t::_prime(def, func, n1, n2);
  }
};

TEST(BicubicSplines, evaluate_cubic_polynom) {
  size_t N = 5;
  auto low = 1.f;
  auto high = 10.f;
  auto def = spline_def_t();
  auto func = [](double x1, double x2) { return x1 * x1 * x2 + x2 * x2; };
  def.f = func;
  def.approx_derivates = true;
  def.axis[0] = std::make_unique<cubic_splines::LinAxis<double>>(low, high, N);
  def.axis[1] = std::make_unique<cubic_splines::ExpAxis<double>>(low, high, N);
  auto spline = cubic_splines::Interpolant<spline_t>(std::move(def), "", "");
  std::uniform_real_distribution<double> dis(low, high);
  std::ofstream file;
  file.open("/tmp/test.txt");
  for (int i = 0; i < 100'000; ++i) {
    auto x = dis(gen);
    /* auto y = dis(gen); */
    auto y = 1.23;
    auto f_inter = spline.evaluate(std::array<double, 2>{x, y});
    auto f_func = func(x, y);
    file << x << ", " << f_func << ", " << f_inter << std::endl;
    /* EXPECT_NEAR(f_func, f_inter, std::max(std::abs(f_func) * 1e-3, 1e-3)); */
  }
};

/* TEST(BicubicSplines, evaluate_absolute_value) { */
/*   size_t N = 30; */
/*   auto low = -10.f; */
/*   auto high = 10.f; */
/*   auto def = spline_def_t(); */
/*   auto func = [](double x1, double x2) { return std::abs(x1) + std::abs(x2); }; */
/*   def.f = func; */
/*   def.axis[0] = std::make_unique<cubic_splines::LinAxis<double>>(low, high, N); */
/*   def.axis[1] = std::make_unique<cubic_splines::LinAxis<double>>(low, high, N); */
/*   auto spline = cubic_splines::Interpolant<spline_t>(std::move(def), "", ""); */
/*   std::uniform_real_distribution<double> dis(low, high); */
/*   for (int i = 0; i < 10'000; ++i) { */
/*     auto x = dis(gen); */
/*     auto y = dis(gen); */
/*     EXPECT_NEAR(func(x, y), spline.evaluate(std::array<double, 2>{x, y}), */
/*                 std::max(std::abs(func(x, y)) * 1e-3, 5e-1)); */
/*   } */
/* } */

/* TEST(BicubicSplines, evaluate_heaviside) { */
/*   size_t N = 30; */
/*   auto low = -0.1f; */
/*   auto high = 0.1f; */
/*   auto def = spline_def_t(); */
/*   auto func = [](double x1, double x2) { return ((x1 > 0) && (x2 > 0)); }; */
/*   def.f = func; */
/*   def.axis[0] = std::make_unique<cubic_splines::LinAxis<double>>(low, high, N); */
/*   def.axis[1] = std::make_unique<cubic_splines::LinAxis<double>>(low, high, N); */
/*   auto spline = cubic_splines::Interpolant<spline_t>(std::move(def), "", ""); */
/*   std::uniform_real_distribution<double> dis(low, high); */
/*   for (int i = 0; i < 10'000; ++i) { */
/*     auto x = dis(gen); */
/*     auto y = dis(gen); */
/*     EXPECT_NEAR(func(x, y), spline.evaluate(std::array<double, 2>{x, y}), 1); */
/*   } */
/* } */

/* TEST(BicubicSplines, evaluate_exp_distributed_nodes) { */
/*   size_t N = 30; */
/*   auto low = 1e-5f; */
/*   auto high = 1e0f; */
/*   auto def = spline_def_t(); */
/*   auto func = [](double x1, double x2) { return std::exp(x1 * x1 + x2 * x2); }; */
/*   def.f = func; */
/*   def.axis[0] = std::make_unique<cubic_splines::ExpAxis<double>>(low, high, N); */
/*   def.axis[1] = std::make_unique<cubic_splines::ExpAxis<double>>(low, high, N); */
/*   auto spline = cubic_splines::Interpolant<spline_t>(std::move(def), "", ""); */
/*   std::uniform_real_distribution<double> dis(low, high); */
/*   for (int i = 0; i < 10'000; ++i) { */
/*     auto x = dis(gen); */
/*     auto y = dis(gen); */
/*     EXPECT_NEAR(func(x, y), spline.evaluate(std::array<double, 2>{x, y}), 1); */
/*   } */
/* } */

/* TEST(BicubicSplines, evaluate_log_func_values_and_axis) { */
/*   size_t N = 10; */
/*   auto low = 1e0f; */
/*   auto high = 1e1f; */
/*   auto def = spline_def_t(); */
/*   auto func = [](double x1, double x2) { return std::pow(x1, 10) * std::pow(x2, 10); };
 */
/*   def.f = func; */
/*   def.f_trafo = std::make_unique<cubic_splines::ExpAxis<double>>(1, 0); */
/*   def.axis[0] = std::make_unique<cubic_splines::ExpAxis<double>>(low, high, N); */
/*   def.axis[1] = std::make_unique<cubic_splines::ExpAxis<double>>(low, high, N); */
/*   auto spline = cubic_splines::Interpolant<spline_t>(std::move(def), "", ""); */
/*   std::uniform_real_distribution<double> dis(low, high); */
/*   for (int i = 0; i < 10; ++i) { */
/*     auto x = dis(gen); */
/*     auto y = dis(gen); */
/*     EXPECT_NEAR(func(x, y), spline.evaluate(std::array<double, 2>{x, y}), */
/*                 std::abs(func(x, y)) * 1e-2); */
/*   } */
/* } */

/* TEST(BicubicSplines, prime) { */
/*   size_t N = 11; */
/*   auto low = 1.e0f; */
/*   auto high = 1.e1f; */
/*   auto def = spline_def_t(); */
/*   auto func = [](double x1, double x2) { return x1 * x1 + x2 * x2 + x1 * x2; }; */
/*   auto df_dx = [](double x1, double x2) { */
/*     return std::array<double, 2>{2 * x1 + x2, 2 * x2 + x1}; */
/*   }; */
/*   def.f = func; */
/*   def.axis[0] = std::make_unique<cubic_splines::LinAxis<double>>(low, high, N); */
/*   def.axis[1] = std::make_unique<cubic_splines::LinAxis<double>>(low, high, N); */
/*   auto splineI = cubic_splines::Interpolant<spline_t>(std::move(def), "", ""); */
/*   std::uniform_real_distribution<double> dis(low, high); */
/*   for (int i = 0; i < 10'000; ++i) { */
/*     auto x = dis(gen); */
/*     auto y = dis(gen); */
/*     EXPECT_NEAR(func(x, y), splineI.evaluate(std::array<double, 2>{x, y}), */
/*                 std::abs(func(x, y)) * 1e-2); */
/*     auto f = df_dx(x, y); */
/*     auto sp = splineI.prime(std::array<double, 2>{x, y}); */
/*     EXPECT_NEAR(f[0], sp[0], std::abs(f[0]) * 1e-2); */
/*     EXPECT_NEAR(f[1], sp[1], std::abs(f[1]) * 1e-2); */
/*   } */
/* } */

/* TEST(BicubicSplines, prime_trafo) { */
/*   size_t N = 11; */
/*   auto low = 1.e0f; */
/*   auto high = 1.e1f; */
/*   auto def = spline_def_t(); */
/*   auto func = [](double x1, double x2) { return x1 * x1 + x2 * x2 + x1 * x2; }; */
/*   auto df_dx = [](double x1, double x2) { */
/*     return std::array<double, 2>{2 * x1 + x2, 2 * x2 + x1}; */
/*   }; */
/*   def.f = func; */
/*   def.f_trafo = std::make_unique<cubic_splines::ExpAxis<double>>(1, 0); */
/*   def.axis[0] = std::make_unique<cubic_splines::ExpAxis<double>>(low, high, N); */
/*   def.axis[1] = std::make_unique<cubic_splines::LinAxis<double>>(low, high, N); */
/*   auto spline = cubic_splines::Interpolant<spline_t>(std::move(def), "", ""); */
/*   std::uniform_real_distribution<double> dis(low, high); */
/*   for (int i = 0; i < 10'000; ++i) { */
/*     auto x = dis(gen); */
/*     auto y = dis(gen); */
/*     EXPECT_NEAR(func(x, y), spline.evaluate(std::array<double, 2>{x, y}), */
/*                 std::abs(func(x, y)) * 1e-2); */
/*     auto f = df_dx(x, y); */
/*     auto sp = spline.prime(std::array<double, 2>{x, y}); */
/*     EXPECT_NEAR(f[0], sp[0], std::abs(f[0]) * 1e-2); */
/*     EXPECT_NEAR(f[1], sp[1], std::abs(f[1]) * 1e-2); */
/*   } */
/* } */

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
