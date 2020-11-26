#pragma once

#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#include <boost/serialization/access.hpp>
#include <cassert>
#include <functional>
#include <memory>
#include <vector>

#include "Axis.h"

/**
 * @brief One dimensional cubic splines class. Tables are build from the lower
 * limit zero with a stepsize from one. If a function has an different
 * definition area it will be  transformed with the Axis transformations
 * specified in the CubicSplines::Definition.
 */
class CubicSplines {
  boost::math::interpolators::cardinal_cubic_b_spline<float> spline;

public:
  /**
   * @brief Initialization of an 1-dim cubic spline with an ordered iterable
   * container stored the values to interpolate to the corresponding axis
   */
  template <typename T1, typename T2>
  CubicSplines(T1 val, T2 lower_lim_derivate, T2 upper_lim_derivate)
      : spline(val.data(), val.size(), 0, 1, lower_lim_derivate,
               upper_lim_derivate) {}

  static constexpr size_t N = 1;

  /**
   * @brief Calculate the function value to the given axis. The interpolated
   * value with no knowledge about the transformation which might has choosen.
   */
  float evaluate(float x) const { return spline(x); };

  float prime(float x) const { return spline.prime(x); };

  float double_prime(float x) const { return spline.double_prime(x); };

  struct Definition;

  class Data;

};

/**
 * @brief Properties of an *1-dim* interpolation object.
 */
struct CubicSplines::Definition {
  /**
   * @brief function to evaluate
   */
  std::function<float(float)> f;

  /**
   * @brief trafo of function values
   */
  std::unique_ptr<Axis> f_trafo = nullptr;

  /**
   * @brief trafo of axis
   */
  std::unique_ptr<Axis> axis;
};

/**
 * @brief Storage class to write and load the interpolation tables from  disk.
 * After reading and writing the object will be destructed.
 */
class CubicSplines::Data {
  std::vector<float> y;
  float lower_lim_derivate;
  float upper_lim_derivate;

  friend class boost::serialization::access;
  template <class Archive> void serialize(Archive &ar, const unsigned int) {
    ar &y;
  }

public:
  Data() = default;

  template <typename T>
  Data(T const &_y, float _lower_lim_derivate, float _upper_lim_derivate)
      : y(_y.begin(), _y.end()), lower_lim_derivate(_lower_lim_derivate),
        upper_lim_derivate(_upper_lim_derivate){};

  CubicSplines build() {
    return CubicSplines(y, lower_lim_derivate, upper_lim_derivate);
  };
};
