#pragma once

#include <cerrno>
#include <fstream>
#include <iostream>
#include <string>
#include <system_error>
#include <tuple>

#include "BicubicSplines.h"

template <typename T1, typename T2 = typename T1::Definition>
class InterpolantBuilder {
  T2 def;

public:
  InterpolantBuilder(T2 &&_def) : def(std::move(_def)){};

  T1 build(){};

  T1 load(std::string path, size_t hash) {
    throw std::system_error(ENOENT, std::generic_category(),
                            "Interpolation tables couldn't be found.");
  };

  template <typename... T> bool save(std::string path, size_t hash, T... args) {
    return false;
  };
};

static std::string SAVE_PATH = "/home/msackel/.local/include/PROPOSAL/tables/";

template <>
BicubicSplines InterpolantBuilder<BicubicSplines>::load(std::string path,
                                                        size_t hash) {
  std::ifstream ifs(SAVE_PATH + "test");
  if (!ifs.is_open()) {
    throw std::system_error(ENOENT, std::generic_category(),
                            "Interpolation tables couldn't be found.");
  }
  std::cout << "load tables" << std::endl;
  auto builder = BicubicSplinesData();
  boost::archive::text_iarchive ia(ifs);
  ia >> builder;
  return builder.build();
};

template <>
template <>
bool InterpolantBuilder<BicubicSplines>::save(std::string path, size_t hash,
                                              Eigen::MatrixXf y,
                                              Eigen::MatrixXf dydx1,
                                              Eigen::MatrixXf dydx2,
                                              Eigen::MatrixXf d2ydx1dx2) {
  auto data = BicubicSplinesData(y, dydx1, dydx2, d2ydx1dx2);
  std::ofstream ofs(SAVE_PATH + "test");
  boost::archive::text_oarchive oa(ofs);
  oa << data;
  return true;
}

template <> BicubicSplines InterpolantBuilder<BicubicSplines>::build() {
  using boost::math::differentiation::finite_difference_derivative;

  try {
    return load(SAVE_PATH, 0);
  } catch (std::system_error const &ex) {
    if (ex.code().value() != ENOENT)
      throw ex;
  }

  auto y = Eigen::MatrixXf(def.nodes[0], def.nodes[1]);
  auto dydx1 = Eigen::MatrixXf(def.nodes[0], def.nodes[1]);
  auto dydx2 = Eigen::MatrixXf(def.nodes[0], def.nodes[1]);
  auto d2ydx1dx2 = Eigen::MatrixXf(def.nodes[0], def.nodes[1]);

  for (size_t n1 = 0; n1 < def.nodes[0]; ++n1) {
    for (size_t n2 = 0; n2 < def.nodes[1]; ++n2) {
      auto x1 = def.x_trafo[0]->back_node(n1);
      auto x2 = def.x_trafo[1]->back_node(n2);

      y(n1, n2) = def.f(x1, x2);
      dydx1(n1, n2) = finite_difference_derivative(
          [this, x2](float x) { return def.f(x, x2); }, x1);
      dydx2(n1, n2) = finite_difference_derivative(
          [this, x1](float x) { return def.f(x1, x); }, x2);
      d2ydx1dx2(n1, n2) = finite_difference_derivative(
          [this, x1, x2](float x_1) {
            return finite_difference_derivative(
                [this, x_1](float x_2) { return def.f(x_1, x_2); }, x2);
          },
          x1);
    }
  }

  bool sucess = save(SAVE_PATH, 0, y, dydx1, dydx2, d2ydx1dx2);
  if (not sucess)
    std::cout << "storage of tables have failed" << std::endl;

  return BicubicSplines(y, dydx1, dydx2, d2ydx1dx2);
}
