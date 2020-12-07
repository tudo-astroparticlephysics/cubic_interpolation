#pragma once

#include <Eigen/Dense>
#include <boost/serialization/access.hpp>
#include <cassert>
#include <memory>
#include <functional>

#include "Axis.h"

namespace cubic_splines {
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

  unsigned int calculate_node(float val, unsigned int n_max) const;

  Eigen::Vector4f build_exponent_vector(float x) const;

  float _evaluate(float x0, float x1) const;

  std::array<float, 2> _prime(float x0, float x1) const;

  float _double_prime(float x0, float x1) const;

public:
  /**
   * @brief Initialization of an 2-dim cubic spline with an Eigen::Xdmatrix
   * stored with the values to the corresponding axis. Beside the function
   * values, the first and second derivates are required.
   */
  BicubicSplines(matrix_ref_t _y, matrix_ref_t _dydx1, matrix_ref_t _dydx2,
                 matrix_ref_t _d2ydx1dx2);

  static constexpr size_t N = 2;

  /**
   * @brief Calculate the function value to the given axis. The interpolated
   * value with no knowledge about the transformation which might has choosen.
   * Requires an iterable container with the *x_i* values stored.
   */
  template <typename T> auto evaluate(T x) const {
    assert(2 == x.size());
    return _evaluate(x[0], x[1]);
  }

  template <typename T> std::array<float, 2> prime(T x) const {
    assert(2 == x.size());
    return _prime(x[0], x[1]);
  }

  template <typename T> float double_prime(T x) const {
    assert(2 == x.size());
    return _double_prime(x[0], x[1]);
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
  std::unique_ptr<cubic_splines::Axis> f_trafo = nullptr;

  /**
   * @brief trafo of axis
   */
  std::array<std::unique_ptr<cubic_splines::Axis>, N> axis;
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
       matrix_ref_t _d2ydx1dx2);

  BicubicSplines build();
};
} // namespace cubic_splines
