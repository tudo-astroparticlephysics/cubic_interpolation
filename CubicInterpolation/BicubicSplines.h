#pragma once

#include <Eigen/Dense>
#include <boost/math/differentiation/finite_difference.hpp>
#include <boost/serialization/access.hpp>
#include <cassert>
#include <functional>

#include "Axis.h"

/**
 * @brief Two dimensional cubic splines class. Tables are build from the lower
 * limit zero with a stepsize from one. If a function has an different
 * definition area it will be  transformed with the Axis transformations
 * specified in the BicubicSplines::Definition.
 */
class BicubicSplines {
  using matrix_t = Eigen::MatrixXf;
  using matrix_ref_t = Eigen::Ref<matrix_t>;

  const Eigen::Matrix4f m, m_t;
  matrix_t y, dydx1, dydx2, d2ydx1dx2;

  template <typename T1, typename T2>
  T2 calculate_node(T1 val, T2 n_max) const {
    auto n = (T2)std::floor(val);
    if (n < 0) {
      return 0;
    } else if (n > n_max - 2u) {
      return n_max - 2u;
    }
    return n;
  }

public:
  /**
   * @brief Initialization of an 2-dim cubic spline with an Eigen::Xdmatrix
   * stored with the values to the corresponding axis. Beside the function
   * values, the first and second derivates are required.
   */
  BicubicSplines(matrix_ref_t _y, matrix_ref_t _dydx1, matrix_ref_t _dydx2,
                 matrix_ref_t _d2ydx1dx2)
      : m((Eigen::Matrix4f() << 1, 0, -3, 2, 0, 0, 3, -2, 0, 1, -2, 1, 0, 0, -1,
           1)
              .finished()),
        m_t(m.transpose()), y(_y), dydx1(_dydx1), dydx2(_dydx2),
        d2ydx1dx2(_d2ydx1dx2) {}

  static constexpr size_t N = 2;

  /**
   * @brief Calculate the function value to the given axis. The interpolated
   * value with no knowledge about the transformation which might has choosen.
   * Requires an iterable container with the *x_i* values stored.
   */
  template <typename T> float evaluate(T x) const {
    assert(2 == x.size());
    auto n0 = calculate_node(x[0], y.cols());
    auto n1 = calculate_node(x[1], y.rows());
    x[0] -= n0;
    x[1] -= n1;
    auto temp = Eigen::Matrix4f(4, 4);
    temp.block<2, 2>(0, 0) = y.block(n0, n1, 2, 2);
    temp.block<2, 2>(2, 0) = dydx1.block(n0, n1, 2, 2);
    temp.block<2, 2>(0, 2) = dydx2.block(n0, n1, 2, 2);
    temp.block<2, 2>(2, 2) = d2ydx1dx2.block(n0, n1, 2, 2);
    auto x1 = Eigen::Vector4f(4);
    auto x2 = Eigen::Vector4f(4);
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

  template <typename T> std::array<float, 2> prime(T x) const {
    using boost::math::differentiation::finite_difference_derivative;
    auto grad = std::array<float, 2>();
    grad[0] = finite_difference_derivative(
        [this, x](float x0) {
          return evaluate(std::array<float, 2>{x0, x[1]});
        },
        x[0]);
    grad[1] = finite_difference_derivative(
        [this, x](float x1) {
          return evaluate(std::array<float, 2>{x[0], x1});
        },
        x[1]);
    return grad;
  }

  template <typename T> float double_prime(T x) const {
    using boost::math::differentiation::finite_difference_derivative;
    return finite_difference_derivative(
        [this, x](float x_1) {
          return finite_difference_derivative(
              [this, x, x_1](float x_2) {
                return evaluate(std::array<float, 2>{x_1, x_2});
              },
              x[1]);
        },
        x[0]);
  }

  struct Definition;

  struct Data;
};

/**
 * @brief Properties of an *2-dim* interpolation object.
 */
struct BicubicSplines::Definition {
  /**
   * @brief function to evaluate
   */
  std::function<float(float, float)> f;

  /**
   * @brief trafo of function values
   */
  std::unique_ptr<Axis> f_trafo = nullptr;

  /**
   * @brief trafo of axis
   */
  std::array<std::unique_ptr<Axis>, N> axis;
};

/**
 * @brief Storage class to write and load the interpolation tables from  disk.
 * After reading and writing the object will be destructed.
 */
class BicubicSplines::Data {

  std::array<size_t, 2> size;
  std::vector<float> y, dydx1, dydx2, d2ydx1dx2;

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
  Data(matrix_ref_t _y, matrix_ref_t _dydx1, matrix_ref_t _dydx2,
       matrix_ref_t _d2ydx1dx2)
      : size({(size_t)_y.rows(), (size_t)_y.cols()}),
        y(std::vector<float>(_y.data(), _y.data() + size[0] * size[1])),
        dydx1(std::vector<float>(_dydx1.data(),
                                 _dydx1.data() + size[0] * size[1])),
        dydx2(std::vector<float>(_dydx2.data(),
                                 _dydx2.data() + size[0] * size[1])),
        d2ydx1dx2(std::vector<float>(_d2ydx1dx2.data(),
                                     _d2ydx1dx2.data() + size[0] * size[1])){};

  BicubicSplines build() {
    return BicubicSplines(
        Eigen::Map<matrix_t>(y.data(), size[0], size[1]),
        Eigen::Map<matrix_t>(dydx1.data(), size[0], size[1]),
        Eigen::Map<matrix_t>(dydx2.data(), size[0], size[1]),
        Eigen::Map<matrix_t>(d2ydx1dx2.data(), size[0], size[1]));
  }
};
