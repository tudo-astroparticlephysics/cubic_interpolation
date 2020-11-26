#pragma once

#include <Eigen/Dense>
#include <cerrno>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <system_error>
#include <tuple>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>

#include <boost/filesystem.hpp>

#include "BicubicSplines.h"
#include "CubicSplines.h"

template <typename T1, typename T2 = typename T1::Definition,
          typename T3 = typename T1::Data>
class InterpolantBuilder {

  T1 load(boost::filesystem::path path, boost::filesystem::path filename) {
    if (boost::filesystem::is_regular_file(path / filename)) {
      std::ifstream ifs((path / filename).c_str());
      auto Data = T3();
      if (ifs.is_open()) {
        boost::archive::text_iarchive ia(ifs);
        ia >> Data;
        return Data.build();
      }
    }
    throw std::system_error(ENOENT, std::generic_category(),
                            "Interpolation tables couldn't be found");
  };

  template <typename... T>
  bool save(boost::filesystem::path path, boost::filesystem::path filename,
            T... args) {
    if (not boost::filesystem::exists(path / filename)) {
      auto data = T3(args...);
      std::ofstream ofs((path / filename).c_str());
      while (ofs.good()) {
        boost::archive::text_oarchive oa(ofs);
        oa << data;
        return true;
      }
    }
    return false;
  };

  template <typename... T>
  auto transform(T2 const& def, T ...args) {
    auto fx = def.f(args...);
    if (def.f_trafo)
      fx = def.f_trafo->transform(fx);
    return fx;
  };

public:
  InterpolantBuilder() = default;

  T1 build(T2 const &def, std::string save_path, std::string filename);
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
  auto x1nodes = def.axis[0]->required_nodes();
  auto x2nodes = def.axis[1]->required_nodes();
  auto y = Eigen::MatrixXf(x1nodes, x2nodes);
  auto dydx1 = Eigen::MatrixXf(x1nodes, x2nodes);
  auto dydx2 = Eigen::MatrixXf(x1nodes, x2nodes);
  auto d2ydx1dx2 = Eigen::MatrixXf(x1nodes, x2nodes);
  auto func = [this, &def](float x1, float x2){return transform(def, x1, x2);};
  for (size_t n1 = 0; n1 < x1nodes; ++n1) {
    for (size_t n2 = 0; n2 < x2nodes; ++n2) {
      auto x1 = def.axis[0]->back_transform(n1);
      auto x2 = def.axis[1]->back_transform(n2);

      auto dfdx1 = def.axis[0]->derive(x1);
      auto dfdx2 = def.axis[1]->derive(x2);

      y(n1, n2) = func(x1, x2);
      dydx1(n1, n2) =
          finite_difference_derivative(
              [this, &func, x2](float x) { return func(x, x2); }, x1) *
          dfdx1;
      dydx2(n1, n2) =
          finite_difference_derivative(
              [this, &func, x1](float x) { return func(x1, x); }, x2) *
          dfdx2;
      d2ydx1dx2(n1, n2) = finite_difference_derivative(
                              [this, &func, x1, x2, dfdx1, dfdx2](float x_1) {
                                return finite_difference_derivative(
                                           [this, &func, x_1, dfdx2](float x_2) {
                                             return func(x_1, x_2);
                                           },
                                           x2) *
                                       dfdx2;
                              },
                              x1) *
                          dfdx1;
    }
  }
  bool sucess = save(save_path, filename, y, dydx1, dydx2, d2ydx1dx2);
  if (not sucess)
    std::cout << "storage of tables have failed" << std::endl;
  return BicubicSplines(y, dydx1, dydx2, d2ydx1dx2);
}

#include <boost/math/differentiation/finite_difference.hpp>

template <>
CubicSplines
InterpolantBuilder<CubicSplines>::build(CubicSplines::Definition const &def,
                                        std::string path,
                                        std::string filename) {
  using boost::math::differentiation::finite_difference_derivative;
  try {
    return load(path, filename);
  } catch (std::system_error const &ex) {
    if (ex.code().value() != ENOENT)
      throw ex;
  }
  auto y = std::vector<float>(def.axis->required_nodes());
  auto func = [this, &def](double x){return transform(def, x);};
  /* auto func = [&def](float x) { */
  /*   auto fx = def.f(x); */
  /*   if (def.f_trafo) */
  /*     fx = def.f_trafo->transform(fx); */
  /*   return fx; */
  /* }; */
  for (size_t n = 0; n < y.size(); ++n)
    y[n] = func(def.axis->back_transform(n));
  auto low = def.axis->back_transform(0);
  auto low_lim_derivate =
      finite_difference_derivative(func, low) * def.axis->derive(low);
  auto up = def.axis->back_transform(y.size() - 1);
  auto up_lim_derivate = finite_difference_derivative(
                             func, def.axis->back_transform(y.size() - 1)) *
                         def.axis->derive(up);
  bool sucess = save(path, filename, y, low_lim_derivate, up_lim_derivate);
  if (not sucess)
    std::cout << "storage of tables have failed" << std::endl;
  return CubicSplines(y, low_lim_derivate, up_lim_derivate);
}
