#include <CubicSplines.h>
#include <Utility.h>
#include <array>
#include <chrono>
#include <iostream>
#include <random>

float func(float x) { return x * x + x; }

int main(int argc, char *argv[]) {

  auto def = CubicSplines::Definition();
  def.f = func;
  def.axis = std::make_unique<LinAxis>(-1.f, 1.f, (size_t)101);

  auto inter = Interpolant<CubicSplines>(
      std::move(def), "/home/msackel/.local/share/PROPOSAL/", "43.txt");

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<float> dis(-1.0, 1.0);

  auto res = 0.f;
  auto point = 0.f;
  auto start = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < 1; ++i) {
    point = dis(gen);

    res = inter.evaluate(point);
  }
  auto stop = std::chrono::high_resolution_clock::now();

  std::cout << point << " * " << point << " + " << point << "= " << res
            << std::endl;

  std::cout << "cubic splien interpolation takes "
            << std::chrono::duration_cast<std::chrono::microseconds>(stop -
                                                                     start)
                   .count()
            << " micro seconds" << std::endl;

  return 0;
}
