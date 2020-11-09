#pragma once

#include <cerrno>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <system_error>
#include <tuple>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "BicubicSplines.h"

template <typename T1, typename T2 = typename T1::Definition,
          typename T3 = typename T1::Data>
class InterpolantBuilder {
  /* T2 def; */

public:
  /* InterpolantBuilder(T2 &&_def) : def(std::move(_def)){}; */
  InterpolantBuilder() = default;

  T1 build(T2 const &def, std::string save_path, std::string filename){};

  T1 load(std::string path, std::string filename) {
    std::ifstream ifs(path + filename);
    if (!ifs.is_open()) {
      ifs.close();
      throw std::system_error(ENOENT, std::generic_category(),
                              "Interpolation tables couldn't be found.");
    }
    auto Data = T3();
    boost::archive::text_iarchive ia(ifs);
    ia >> Data;
    return Data.build();
  };

  template <typename... T>
  bool save(std::string path, std::string filename, T... args) {
    auto data = T3(args...);
    std::ofstream ofs(path + filename);
    while (ofs.good()) {
      boost::archive::text_oarchive oa(ofs);
      oa << data;
      return true;
    }
    return false;
  };
};

template <>
BicubicSplines
InterpolantBuilder<BicubicSplines>::build(BicubicSplines::Definition const &def,
                                          std::string save_path,
                                          std::string filename) {
  using boost::math::differentiation::finite_difference_derivative;

  try {
    return load(save_path, filename);
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
          [this, &def, x2](float x) { return def.f(x, x2); }, x1);
      dydx2(n1, n2) = finite_difference_derivative(
          [this, &def, x1](float x) { return def.f(x1, x); }, x2);
      d2ydx1dx2(n1, n2) = finite_difference_derivative(
          [this, &def, x1, x2](float x_1) {
            return finite_difference_derivative(
                [this, &def, x_1](float x_2) { return def.f(x_1, x_2); }, x2);
          },
          x1);
    }
  }

  bool sucess = save(save_path, filename, y, dydx1, dydx2, d2ydx1dx2);
  if (not sucess)
    std::cout << "storage of tables have failed" << std::endl;

  return BicubicSplines(y, dydx1, dydx2, d2ydx1dx2);
}
