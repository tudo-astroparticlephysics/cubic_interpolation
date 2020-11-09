#include <BicubicSplines.h>
#include <Utility.h>
#include <array>
#include <chrono>
#include <iostream>
#include <random>

float func(float x_1, float x_2) { return x_1 * x_1 + x_2 * x_2 + x_1 * x_2; }

int main(int argc, char *argv[]) {
  constexpr size_t size = 11;

  auto low = std::array<float, 2>{-1, -1};
  auto high = std::array<float, 2>{1, 1};
  auto step_size = std::array<float, 2>{(high[0] - low[0]) / (size - 1),
                                        (high[0] - low[0]) / (size - 1)};

  auto def = BicubicSplines::Definition();
  def.f = func;
  def.axis[0] = std::make_unique<LinAxis>(low[0], high[0], size);
  def.axis[1] = std::make_unique<LinAxis>(low[1], high[1], size);

  auto inter = Interpolant<BicubicSplines>(
      std::move(def), "/home/msackel/.local/share/PROPOSAL/", "42.txt");

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<float> dis(-1.0, 1.0);

  auto res = 0.f;
  auto point = std::array<float, 2>();
  auto start = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < 1'000; ++i) {
    point[0] = dis(gen);
    point[1] = dis(gen);

    res = inter.evaluate(point);
  }
  auto stop = std::chrono::high_resolution_clock::now();

  std::cout << point[0] << " * " << point[0] << " + " << point[1] << " * "
            << point[1] << " + " << point[0] << " * " << point[1] << " = "
            << res << std::endl;

  std::cout << "cubic splien interpolation takes "
            << std::chrono::duration_cast<std::chrono::microseconds>(stop -
                                                                     start)
                   .count()
            << " micro seconds" << std::endl;

  return 0;
}
