#pragma once

#include <Eigen/Dense>
#include <boost/math/differentiation/finite_difference.hpp>
#include <boost/serialization/access.hpp>
#include <cassert>
#include <functional>

#include "Axis.h"

class BicubicSplines {

  const Eigen::Matrix4f m;
  const Eigen::Matrix4f m_t;

  using matrix_t = Eigen::MatrixXf;
  matrix_t y;
  matrix_t dydx1;
  matrix_t dydx2;
  matrix_t d2ydx1dx2;

public:
  BicubicSplines() = default;
  BicubicSplines(Eigen::Ref<matrix_t>, Eigen::Ref<matrix_t>,
                 Eigen::Ref<matrix_t>, Eigen::Ref<matrix_t>);

  static constexpr size_t N = 2;

  struct Definition;

  struct Data;

  template <typename T, typename Y> float evaluate(T node, Y x);
};

struct BicubicSplines::Definition {
  std::function<float(float, float)> f;
  std::unique_ptr<Axis> f_trafo = nullptr;
  std::array<std::unique_ptr<Axis>, N> axis;
};

class BicubicSplines::Data {

  std::array<size_t, 2> size;
  std::vector<float> y;
  std::vector<float> dydx1;
  std::vector<float> dydx2;
  std::vector<float> d2ydx1dx2;

  friend class boost::serialization::access;
  template <class Archive> void serialize(Archive &ar, const unsigned int) {
    ar &size;
    ar &y;
    ar &dydx1;
    ar &dydx2;
    ar &d2ydx1dx2;
  }

public:
  Data() = default;
  Data(matrix_t _y, matrix_t _dydx1, matrix_t _dydx2, matrix_t _d2ydx1dx2);

  BicubicSplines build();
};

BicubicSplines::BicubicSplines(Eigen::Ref<matrix_t> _y,
                               Eigen::Ref<matrix_t> _dydx1,
                               Eigen::Ref<matrix_t> _dydx2,
                               Eigen::Ref<matrix_t> _d2ydx1dx2)
    : m((Eigen::Matrix4f() << 1, 0, -3, 2, 0, 0, 3, -2, 0, 1, -2, 1, 0, 0, -1,
         1)
            .finished()),
      m_t(m.transpose()), y(_y), dydx1(_dydx1), dydx2(_dydx2),
      d2ydx1dx2(_d2ydx1dx2) {}

template <typename T1, typename T2>
float BicubicSplines::evaluate(T1 node, T2 x) {
  assert(2 == node.size());
  assert(2 == x.size());

  Eigen::Matrix4f temp(4, 4);
  temp.block<2, 2>(0, 0) = y.block(node[0], node[1], 2, 2);
  temp.block<2, 2>(2, 0) = dydx1.block(node[0], node[1], 2, 2);
  temp.block<2, 2>(0, 2) = dydx2.block(node[0], node[1], 2, 2);
  temp.block<2, 2>(2, 2) = d2ydx1dx2.block(node[0], node[1], 2, 2);

  Eigen::Vector4f x1(4);
  Eigen::Vector4f x2(4);
  for (size_t i = 0; i < 4; ++i) {
    x1(i) = 1;
    x2(i) = 1;
    for (size_t j = 1; j <= i; ++j) {
      x1(i) *= x[0];
      x2(i) *= x[1];
    }
  }

  return x1.dot((m_t * (temp * m)) * x2);
}

BicubicSplines::Data::Data(matrix_t _y, matrix_t _dydx1, matrix_t _dydx2,
                           matrix_t _d2ydx1dx2)
    : size({(size_t)_y.rows(), (size_t)_y.cols()}),
      y(std::vector<float>(_y.data(), _y.data() + size[0] * size[1])),
      dydx1(
          std::vector<float>(_dydx1.data(), _dydx1.data() + size[0] * size[1])),
      dydx2(
          std::vector<float>(_dydx2.data(), _dydx2.data() + size[0] * size[1])),
      d2ydx1dx2(std::vector<float>(_d2ydx1dx2.data(),
                                   _d2ydx1dx2.data() + size[0] * size[1])){};

BicubicSplines BicubicSplines::Data::build() {
  using Eigen::Map;
  return BicubicSplines(Map<matrix_t>(y.data(), size[0], size[1]),
                        Map<matrix_t>(dydx1.data(), size[0], size[1]),
                        Map<matrix_t>(dydx2.data(), size[0], size[1]),
                        Map<matrix_t>(d2ydx1dx2.data(), size[0], size[1]));
};
