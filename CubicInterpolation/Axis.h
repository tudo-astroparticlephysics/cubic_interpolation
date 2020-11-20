/**
 * @file dlfajsldfjkasl;d dsfkld
 *
 */

#pragma once

#include <cassert>
#include <cmath>
#include <iostream>
#include <ostream>

/**
 * @brief Limits the definition area of the function which will be stored.
 * Choose the axis transformation wisely to keep the number of support points
 * and the approximation error low. Evaluate points out of range will cause
 * undefined behavour.
 */
class Axis {
protected:
  float low, high, stepsize;

  virtual void print(std::ostream &os) const {
    os << "low: " << low << ", high: " << high << ", stepsize: " << stepsize;
  };

  friend std::ostream &operator<<(std::ostream &out, const Axis &);

public:
  /**
   * @brief Initialize the minimum required variables to build an Interpolation
   * Axis.
   *
   * @param _low lower limit of the axis
   * @param _high upper limit of the axis
   * @param _stepsize stepsize in values of the transformation
   */
  Axis(float _low, float _high, float _stepsize)
      : low(_low), high(_high), stepsize(_stepsize) {}

  /**
   * @brief Return the corresponding lower node which is nearest to the point
   * and manipulate the *x*-value to the relative distance from the node to the
   * next one.
   */
  virtual float transform(float x) const noexcept = 0;

  /**
   * @brief Calculate the corresponding Point to the *n*-th node. The *0*-th
   * node correspond to the lower limit. The *N-1*-th represent the upper limit
   * of the axis.
   */
  virtual float back_transform(float x) const noexcept = 0;

  /**
   * @brief Calculates the number of nodes required for generating the
   * interpolation table.
   */
  virtual size_t required_nodes() const {
    auto nodes = transform(high);
    if (std::floor(nodes) == nodes)
      return nodes + 1;
    return nodes + 2;
  }

  virtual float derive(float x) const = 0;
};

/**
 * @brief Linear axis to describe data which varify not to much in orders of
 * magnitudes. It's the fastest axis evaluation.
 */
class LinAxis : public Axis {
  void print(std::ostream &os) const {
    os << "LinAxis(";
    Axis::print(os);
    os << ")";
  }

public:
  /**
   * @brief Linear Axis initalized with stepsize.
   */
  LinAxis(float _low, float _high, float _stepsize)
      : Axis(_low, _high, _stepsize) {}

  /**
   * @brief Linear Axis with number of nodes.
   */
  LinAxis(float _low, float _high, size_t _nodes) : Axis(_low, _high, 0.f) {
    stepsize = (high - low) / (_nodes - 1);
  }

  inline float transform(float x) const noexcept final {
    return (x - low) / stepsize;
  }
  inline float back_transform(float x) const noexcept final {
    return x * stepsize + low;
  }

  float derive(float) const final { return stepsize; }
};

/**
 * @brief Exponential axis to interpolate over many different scales. As basis
 * the euler number choosen by default if no further specified. Nodes will
 * distributed in form of \f$ \exp(n \cdot \text{stepsize}) \f$
 */
class ExpAxis : public Axis {
  void print(std::ostream &os) const {
    os << "ExpAxis(";
    Axis::print(os);
    os << ")";
  }

public:
  /**
   * @brief Exponential Axis initialized with stepsize.
   */
  ExpAxis(float _low, float _high, float _stepsize = 1)
      : Axis(_low, _high, _stepsize) {}

  /**
   * @brief Exponential Axis initialized with number of nodes.
   */
  ExpAxis(float _low, float _high, size_t _nodes) : Axis(_low, _high, 0.f) {
    stepsize = std::log(high / low) / (_nodes - 1);
  }

  inline float transform(float x) const noexcept final {
    return std::log(x / low) / stepsize;
  }

  inline float back_transform(float x) const noexcept final {
    return low * std::exp(x * stepsize);
  }

  float derive(float x) const final { return (x * stepsize); }
};

std::ostream &operator<<(std::ostream &out, const Axis &axis) {
  axis.print(out);
  return out;
}
