#include <cstdlib>
#include <iostream>
#include <memory>

#include "CubicSplines.h"
#include "Utility.h"

float func(float x) { return x * x + x; }

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cout << "Usage: 1dim_example POINT\n"
              << "\n"
              << "  POINT     point to evaluate the function\n"
              << std::endl;
    return 0;
  }

  auto def = CubicSplines::Definition();
  def.f = func;
  def.axis = std::make_unique<LinAxis>(-4.f, 4.f, (size_t)20);

  auto path = "/home/msackel/.local/share/PROPOSAL/";
  auto tablename = "43.txt";
  auto inter = Interpolant<CubicSplines>(std::move(def), path, tablename);

  auto point = std::atof(argv[1]);
  auto res = inter.evaluate(point);

  std::cout << "f(" << point << "): " << res << std::endl;

  return 0;
}
