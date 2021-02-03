#pragma once

#include "Axis.h"

#include <array>
#include <functional>
#include <memory>
#include <vector>

namespace cubic_splines {

/**
 * @brief Two dimensional cubic splines class. Tables are build from the lower
 * limit zero with a stepsize from one. If a function has an different
 * definition area it will be  transformed with the Axis transformations
 * specified in the BicubicSplines::Definition.
 */
template <typename T> class BicubicSplines {
public:
  using type = T;
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
    std::function<T(T, T)> f;                     // function to evaluate
    std::unique_ptr<Axis<T>> f_trafo;             // trafo of function values
    std::array<std::unique_ptr<Axis<T>>, N> axis; // trafo of axis
    bool approx_derivates;

    const std::array<std::unique_ptr<Axis<T>>, N> &GetAxis() const { return axis; };
  };

  BicubicSplines(Definition const &);
  BicubicSplines(Definition const &, std::string, std::string);

protected:
  BicubicSplines(RuntimeData);

  std::shared_ptr<RuntimeData> data;

  std::tuple<T, T> back_transform(Definition const &, long unsigned int,
                                  long unsigned int) const;

  template <typename T1>
  std::array<T, 2> _prime(T1 func, std::array<std::unique_ptr<Axis<T>>, 2> const &axis,
                          unsigned int n1, unsigned int n2);

  /* template <typename T1> */
  /* std::vector<std::tuple<unsigned int, T>> _prime(std::vector<T> const &yi, T1 func, */
  /*                                                 unsigned int n_max); */

  template <typename T1>
  T _double_prime(Definition const &def, T1 func, unsigned int n1, unsigned int n2);

public:
  /**
   * @brief Calculate the function value to the given axis. The interpolated
   * value with no knowledge about the transformation which might has choosen.
   * Requires an iterable container with the *x_i* values stored.
   */
  T evaluate(T x0, T x1) const;

  template <typename T1> auto evaluate(T1 iterable) const {
    return evaluate(iterable[0], iterable[1]);
  }

  std::array<T, 2> prime(T x0, T x1) const;

  template <typename T1> auto prime(T1 iterable) const {
    return prime(iterable[0], iterable[1]);
  }

  T double_prime(T x0, T x1) const;

  template <typename T1> auto double_prime(T1 iterable) const {
    return double_prime(iterable[0], iterable[1]);
  }
};
} // namespace cubic_splines
