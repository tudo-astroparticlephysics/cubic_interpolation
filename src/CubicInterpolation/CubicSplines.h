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
template <typename T> class CubicSplines {
public:
  using type = T;

  struct StorageData;

  struct RuntimeData;

  static constexpr size_t N = 1;

  /**
   * @brief Properties of an *1-dim* interpolation object.
   */
  struct Definition {
    std::function<T(T)> f;                                     // function to evaluate
    std::unique_ptr<cubic_splines::Axis<T>> f_trafo = nullptr; // trafo of function values
    std::unique_ptr<cubic_splines::Axis<T>> axis;              // trafo of axis

    const Axis<T> &GetAxis() const { return *axis; };
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
  T evaluate(T x) const;

  T prime(T x) const;

  T double_prime(T x) const;
};


} // namespace cubic_splines
