#pragma once

#include <functional>
#include <memory>

#include "Axis.h"

namespace cubic_splines {
/**
 * @brief One dimensional cubic splines class. Tables are build from the lower
 * limit zero with a stepsize from one. If a function has an different
 * definition area it will be  transformed with the Axis transformations
 * specified in the CubicSplines::Definition.
 */
class CubicSplines {

public:
  struct StorageData;

  struct RuntimeData;

  static constexpr size_t N = 1;

  /**
   * @brief Properties of an *1-dim* interpolation object.
   */
  struct Definition {
    std::function<float(float)> f;                          // function to evaluate
    std::unique_ptr<cubic_splines::Axis> f_trafo = nullptr; // trafo of function values
    std::unique_ptr<cubic_splines::Axis> axis;              // trafo of axis

    const Axis &GetAxis() const { return *axis; };
  };

  CubicSplines(Definition const &);
  CubicSplines(Definition const &, std::string, std::string);

private:
  CubicSplines(RuntimeData);

  std::shared_ptr<RuntimeData> data;

public:
  /**
   * @brief Calculate the function value to the given axis. The interpolated
   * value with no knowledge about the transformation which might has choosen.
   */
  float evaluate(float x) const;

  float prime(float x) const;

  float double_prime(float x) const;
};
} // namespace cubic_splines
