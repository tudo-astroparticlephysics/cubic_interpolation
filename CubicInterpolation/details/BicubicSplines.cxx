#include "../BicubicSplines.h"

#include <boost/math/differentiation/finite_difference.hpp>

namespace cubic_splines {

unsigned int BicubicSplines::calculate_node(float val,
                                            unsigned int n_max) const {
  auto n = std::floor(val);
  if (n < 0) {
    return 0;
  } else if (n > n_max - 2u) {
    return n_max - 2u;
  }
  return n;
}

Eigen::Vector4f BicubicSplines::build_exponent_vector(float x) const {
  auto x_vec = Eigen::Vector4f(4);
  for (size_t i = 0; i < 4; ++i) {
    x_vec(i) = 1;
    for (size_t j = 1; j <= i; ++j)
      x_vec(i) *= x;
  }
  return x_vec;
}

float BicubicSplines::_evaluate(float x0, float x1) const {
  auto n0 = calculate_node(x0, y.cols());
  auto n1 = calculate_node(x1, y.rows());

  auto temp = Eigen::Matrix4f(4, 4);
  temp.block<2, 2>(0, 0) = y.block(n0, n1, 2, 2);
  temp.block<2, 2>(2, 0) = dydx1.block(n0, n1, 2, 2);
  temp.block<2, 2>(0, 2) = dydx2.block(n0, n1, 2, 2);
  temp.block<2, 2>(2, 2) = d2ydx1dx2.block(n0, n1, 2, 2);

  auto v1 = build_exponent_vector(x0 - n0);
  auto v2 = build_exponent_vector(x1 - n1);

  return v1.dot((m_t * (temp * m)) * v2);
}

std::array<float, 2> BicubicSplines::_prime(float x0, float x1) const {
  using boost::math::differentiation::finite_difference_derivative;
  auto grad = std::array<float, 2>();
  grad[0] = finite_difference_derivative(
      [this, x1](float x_0) {
        return evaluate(std::array<float, 2>{x_0, x1});
      },
      x0);
  grad[1] = finite_difference_derivative(
      [this, x0](float x_1) {
        return evaluate(std::array<float, 2>{x0, x_1});
      },
      x1);
  return grad;
}

float BicubicSplines::_double_prime(float x0, float x1) const {
  using boost::math::differentiation::finite_difference_derivative;
  return finite_difference_derivative(
      [this, x0, x1](float x_1) {
        return finite_difference_derivative(
            [this, x1, x_1](float x_2) {
              return evaluate(std::array<float, 2>{x_1, x_2});
            },
            x1);
      },
      x0);
}

BicubicSplines::BicubicSplines(matrix_ref_t _y, matrix_ref_t _dydx1,
                               matrix_ref_t _dydx2, matrix_ref_t _d2ydx1dx2)
    : m((Eigen::Matrix4f() << 1, 0, -3, 2, 0, 0, 3, -2, 0, 1, -2, 1, 0, 0, -1,
         1)
            .finished()),
      m_t(m.transpose()), y(_y), dydx1(_dydx1), dydx2(_dydx2),
      d2ydx1dx2(_d2ydx1dx2) {}

BicubicSplines::Data::Data(matrix_ref_t _y, matrix_ref_t _dydx1,
                           matrix_ref_t _dydx2, matrix_ref_t _d2ydx1dx2)
    : size({(size_t)_y.rows(), (size_t)_y.cols()}),
      y(std::vector<float>(_y.data(), _y.data() + size[0] * size[1])),
      dydx1(
          std::vector<float>(_dydx1.data(), _dydx1.data() + size[0] * size[1])),
      dydx2(
          std::vector<float>(_dydx2.data(), _dydx2.data() + size[0] * size[1])),
      d2ydx1dx2(std::vector<float>(_d2ydx1dx2.data(),
                                   _d2ydx1dx2.data() + size[0] * size[1])){};

BicubicSplines BicubicSplines::Data::build() {
  return BicubicSplines(
      Eigen::Map<matrix_t>(y.data(), size[0], size[1]),
      Eigen::Map<matrix_t>(dydx1.data(), size[0], size[1]),
      Eigen::Map<matrix_t>(dydx2.data(), size[0], size[1]),
      Eigen::Map<matrix_t>(d2ydx1dx2.data(), size[0], size[1]));
}
} // namespace cubic_splines
