#pragma once

#include <ostream>

namespace cubic_splines {
/**
 * @brief Limits the definition area of the function which will be stored and
 * praameterize the node spaceing.  Choose the axis transformation wisely to
 * keep the number of support points and the approximation error low. Evaluate
 * points out of range will cause undefined behavour.
 */
class Axis {
protected:
  float low, high, stepsize;

  virtual void print(std::ostream &os) const;

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
  Axis(float _low, float _high, float _stepsize);

  /**
   * @brief Return the corresponding lower node which is nearest to the
   * point and manipulate the *x*-value to the relative distance from the
   * node to the next one.
   */
  virtual float transform(float x) const noexcept = 0;

  /**
   * @brief Calculate the corresponding point to the *n*-th node. The *0*-th
   * node correspond to the lower limit. The *N-1*-th represent the upper limit
   * of the axis. The decimal digits specify the relative distance to the next
   * node.
   */
  virtual float back_transform(float x) const noexcept = 0;

  /**
   * @brief Calculates the required number of nodes.
   */
  virtual size_t required_nodes() const;

  /**
   * @brief Lower axis limit.
   */
  float GetLow() const noexcept { return low; }

  /**
   * @brief Upper axis limit.
   */
  float GetHigh() const noexcept { return high; }

  /**
   * @brief Derivate of the Axis forward transformation.
   */
  virtual float derive(float x) const = 0;
};

std::ostream &operator<<(std::ostream &out, const Axis &axis);
} // namespace cubic_splines

namespace cubic_splines {
/**
 * @brief Exponential axis to interpolate over many different scales. As basis
 * the euler number choosen by default if no further specified. Nodes will
 * distributed in form of \f$ \exp(n \cdot \text{stepsize}) \f$
 */
class ExpAxis : public Axis {

  void print(std::ostream &os) const;

public:
  /**
   * @brief Exponential Axis initialized with stepsize.
   */
  ExpAxis(float _low, float _high, float _stepsize = 1);

  /**
   * @brief Exponential Axis initialized with number of nodes.
   */
  ExpAxis(float _low, float _high, size_t _nodes);

  float transform(float x) const noexcept final;
  float back_transform(float x) const noexcept final;
  float derive(float x) const final;
};
} // namespace cubic_splines

namespace cubic_splines {
/**
 * @brief Linear axis to describe data which varify not to much in orders of
 * magnitudes. It's the fastest axis evaluation.
 */
class LinAxis : public Axis {

  void print(std::ostream &os) const;

public:
  /**
   * @brief Linear Axis initalized with stepsize.
   */
  LinAxis(float _low, float _high, float _stepsize);

  /**
   * @brief Linear Axis with number of nodes.
   */
  LinAxis(float _low, float _high, size_t _nodes);

  float transform(float x) const noexcept final;
  float back_transform(float x) const noexcept final;
  float derive(float) const final;
};
} // namespace cubic_splines
