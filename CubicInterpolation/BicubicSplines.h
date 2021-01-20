#pragma once

#include "Axis.h"

#include <functional>
#include <memory>

namespace cubic_splines {

/**
 * @brief Two dimensional cubic splines class. Tables are build from the lower
 * limit zero with a stepsize from one. If a function has an different
 * definition area it will be  transformed with the Axis transformations
 * specified in the BicubicSplines::Definition.
 */
class BicubicSplines {
public:
  /**
   * @brief Storage class to write and load the interpolation tables from  disk.
   * After reading and writing the object will be destructed.
   */
  struct StorageData;

  struct RuntimeData;

  static constexpr size_t N = 2;

  /**
   * @brief Properties of an *2-dim* interpolation object.
   */
  struct Definition {
    std::function<double(double, double)> f;      // function to evaluate
    std::unique_ptr<Axis> f_trafo = nullptr;   // trafo of function values
    std::array<std::unique_ptr<Axis>, N> axis; // trafo of axis

    const std::array<std::unique_ptr<Axis>, N> &GetAxis() const { return axis; };
  };

  BicubicSplines(Definition const &);
  BicubicSplines(Definition const &, std::string, std::string);

private:
  BicubicSplines(RuntimeData);

  std::shared_ptr<RuntimeData> data;

public:
  /**
   * @brief Calculate the function value to the given axis. The interpolated
   * value with no knowledge about the transformation which might has choosen.
   * Requires an iterable container with the *x_i* values stored.
   */
  double evaluate(double x0, double x1) const;

  template <typename T> auto evaluate(T iterable) const {
    return evaluate(iterable[0], iterable[1]);
  }

  std::array<double, 2> prime(double x0, double x1) const;

  template <typename T> auto prime(T iterable) const {
    return prime(iterable[0], iterable[1]);
  }

  double double_prime(double x0, double x1) const;

  template <typename T> auto double_prime(T iterable) const {
    return double_prime(iterable[0], iterable[1]);
  }
};
} // namespace cubic_splines
