/**
 * @file dlfajsldfjkasl;d dsfkld
 *
 */

#pragma once

#include <cmath>
#include <ostream>

/**
 * @brief Limits the definition area of the function which will be stored.
 * Choose the axis transformation wisely to keep the number of support points
 * and the approximation error low. Evaluate points out of range will cause
 * undefined behavour.
 */
class Axis {
protected:
  float low;
  float high;
  float stepsize;

  virtual void print(std::ostream &) const = 0;

  virtual float transform(float x) const noexcept = 0;
  virtual float back_transform(float x) const noexcept = 0;

  friend std::ostream &operator<<(std::ostream &out, const Axis &);

public:
  /**
   * @brief Initialize the minimum required variables to build an Interpolation
   * Axis.
   *
   * @param _low lower limit of the interpolation tabel
   * @param _high upper limit of the interpolation tabel
   * @param _stepsize stepsize in values of the transformation
   */
  Axis(float _low, float _high, float _stepsize)
      : low(_low), high(_high), stepsize(_stepsize) {}

  /**
   * @brief Calculates the number of nodes required for generating the
   * interpolation table.
   */
  virtual size_t required_nodes() const {
    return std::ceil((high - low) / stepsize) + 1;
  }

  /**
   * @brief Calculate the corresponding lower node which is nearest to the point
   */
  virtual size_t node(float &x) {
    size_t n = std::floor((x - low) / stepsize);
    x -= low + n * stepsize;
    x *= 1 / stepsize;
    return n;
  }

  /**
   * @brief Calculate the corresponding Point to the n-th node.
   */
  virtual float back_node(size_t n) { return low + n * stepsize; }
};

/**
 * @brief Logarithm axis to interpolate over many different scales. As basis is
 * the euler number choosen.
 */
struct LogAxis : public Axis {
  LogAxis(float _low, float _high, float _stepsize)
      : Axis(_low, _high, _stepsize) {}

  LogAxis(float _low, float _high, size_t _nodes) : Axis(_low, _high) {
    transform(high);
    transform(low);
    stepsize = (high - low) / (_nodes - 1);
  }

  inline float transform(float x) const noexcept final { return std::log(x); }

  inline float back_transform(float x) const noexcept final {
    return std::exp(x);
  }

protected:
  void print(std::ostream &) const;
};

struct LinAxis : public Axis {
  LinAxis(float _low, float _high, size_t _nodes) : Axis(_low, _high, 0.f) {
    stepsize = (high - low) / (_nodes - 1);
  }
  LinAxis(float _low, float _high, float _stepsize)
      : Axis(_low, _high, _stepsize) {}
  inline float transform(float x) const noexcept final { return x; }
  inline float back_transform(float x) const noexcept final { return x; }

protected:
  void print(std::ostream &) const;
};

std::ostream &operator<<(std::ostream &out, const Axis &axis) {
  axis.print(out);
  return out;
}

void LogAxis::print(std::ostream &os) const {
  os << "LogAxis(low: " << low << ", high: " << high
     << ", stepsize: " << stepsize << ")";
}

void LinAxis::print(std::ostream &os) const {
  os << "LinAxis(low: " << low << ", high: " << high
     << ", stepsize: " << stepsize << ")";
}
