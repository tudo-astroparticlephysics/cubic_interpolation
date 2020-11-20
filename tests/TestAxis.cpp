#include "Axis.h"
#include "gtest/gtest.h"
#include <cmath>

TEST(LinAxis, transform) {
  auto MIN_VALUE = 0.f;
  auto MAX_VALUE = 10.f;
  auto STEPSIZE = 1.f;
  auto lin = LinAxis(MIN_VALUE, MAX_VALUE, STEPSIZE);

  auto x = 1.f;
  auto x_transformed = lin.transform(x);
  EXPECT_EQ(1, std::floor(x_transformed));
  EXPECT_NEAR(x, x_transformed, 1e-6);

  x = 0.f;
  x_transformed = lin.transform(x);
  EXPECT_EQ(0, std::floor(x_transformed));
  EXPECT_NEAR(x, x_transformed, 1e-6);

  x = 3.45f;
  x_transformed = lin.transform(x);
  EXPECT_EQ(3, std::floor(x_transformed));
  EXPECT_NEAR(x, x_transformed, 1e-6);
}

TEST(LinAxis, back_transform) {
  auto MIN_VALUE = 0.f;
  auto MAX_VALUE = 10.f;
  auto STEPSIZE = 1.f;
  auto lin = LinAxis(MIN_VALUE, MAX_VALUE, STEPSIZE);

  auto n = 3;
  auto x = lin.back_transform(n);
  EXPECT_NEAR(3.f, x, 1e-6);
}

TEST(LinAxis, number_of_nodes) {
  auto MIN_VALUE = 0.f;
  auto MAX_VALUE = 10.f;

  auto stepsize = 0.5f;
  auto lin = LinAxis(MIN_VALUE, MAX_VALUE, stepsize);
  EXPECT_EQ(21, lin.required_nodes());

  stepsize = 0.333f;
  lin = LinAxis(MIN_VALUE, MAX_VALUE, stepsize);
  EXPECT_EQ(32, lin.required_nodes());
}

TEST(ExpAxis, node) {
  auto MIN_VALUE = 1e0f;
  auto MAX_VALUE = 1e14f;
  auto STEPSIZE = std::log(10.f);
  auto log = ExpAxis(MIN_VALUE, MAX_VALUE, STEPSIZE);

  auto x = 1.234e4f;
  auto x_transformed = log.transform(x);
  auto rel_dist = x_transformed - std::floor(x_transformed);
  EXPECT_EQ(4, std::floor(x_transformed));
  EXPECT_NEAR(rel_dist, (std::log(1.234e4f) - std::log(1e4)) / STEPSIZE, 1e-6);
}

TEST(ExpAxis, back_transform) {
  auto MIN_VALUE = 1e0f;
  auto MAX_VALUE = 1e14f;
  auto STEPSIZE = std::log(10.f);
  auto log = ExpAxis(MIN_VALUE, MAX_VALUE, STEPSIZE);

  auto n = 3;
  auto x = log.back_transform(n);
  EXPECT_NEAR(1e3, x, 1e3 * 1e-6);
}

TEST(ExpAxis, number_of_nodes) {
  auto log = ExpAxis(1e0f, 1e14f, std::log(10.f));
  EXPECT_EQ(15, log.required_nodes());

  log = ExpAxis(2, 2048, std::log(2.f));
  EXPECT_EQ(11, log.required_nodes());

  log = ExpAxis(2, 2048, std::log(1.9f));
  EXPECT_EQ(12, log.required_nodes());
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
