#pragma once

#include <ostream>

namespace cubic_splines {
/**
 * @brief Limits the definition area of the function which will be stored and
 * praameterize the node spaceing.  Choose the axis transformation wisely to
 * keep the number of support points and the approximation error low. Evaluate
 * points out of range will cause undefined behavour.
 */
template <typename T> class Axis {
protected:
  T low, high, stepsize;

  virtual void print(std::ostream &os) const {
    os << "low: " << low << ", high: " << high << ", stepsize: " << stepsize;
  };

  /* friend std::ostream &operator<<(std::ostream &out, const Axis<T> &); */

public:
  /**
   * @brief Initialize the minimum required variables to build an Interpolation
   * Axis.
   *
   * @param _low lower limit of the axis
   * @param _high upper limit of the axis
   * @param _stepsize stepsize in values of the transformation
   */
  Axis(T _low, T _high, T _stepsize);
  virtual ~Axis() = default;

  /**
   * @brief Return the corresponding lower node which is nearest to the
   * point and manipulate the *x*-value to the relative distance from the
   * node to the next one.
   */
  virtual T transform(T) const = 0;

  /**
   * @brief Calculate the corresponding point to the *n*-th node. The *0*-th
   * node correspond to the lower limit. The *N-1*-th represent the upper limit
   * of the axis. The decimal digits specify the relative distance to the next
   * node.
   */
  virtual T back_transform(T) const = 0;

  /**
   * @brief Calculates the required number of nodes.
   */
  unsigned int required_nodes() const;

  /**
   * @brief Lower axis limit.
   */
  auto GetLow() const noexcept { return low; }

  /**
   * @brief Upper axis limit.
   */
  auto GetHigh() const noexcept { return high; }

  /**
   * @brief Stepsize of axis.
   */
  auto GetStepsize() const noexcept { return stepsize; }

  /**
   * @brief Derivate of the Axis forward transformation.
   */
  virtual T derive(T x) const = 0;
  virtual T back_derive(T x) const = 0;
};

template <typename T> std::ostream &operator<<(std::ostream &out, const Axis<T> &axis) {
  axis.print(out);
  return out;
}

/**
 * @brief Exponential axis to interpolate over many different scales. As basis
 * the euler number choosen by default if no further specified. Nodes will
 * distributed in form of \f$ \exp(n \cdot \text{stepsize}) \f$
 */
template <typename T> class ExpAxis : public Axis<T> {

  void print(std::ostream &os) const {
    os << "ExpAxis(";
    Axis<T>::print(os);
    os << ")";
  }

public:
  /**
   * @brief Exponential Axis initialized with stepsize.
   */
  ExpAxis(T _low, T _high, T _stepsize = 1);

  /**
   * @brief Exponential Axis initialized with number of nodes.
   */
  ExpAxis(T _low, T _high, size_t _nodes);

  T transform(T x) const final;
  T back_transform(T t) const final;

  T derive(T x) const final;
  T back_derive(T t) const final;
};

/**
 * @brief Exponential axis to interpolate over many different scales. As basis
 * the euler number choosen by default if no further specified. Nodes will
 * distributed in form of \f$ \text{low} \cdot (\exp(n \cdot \text{stepsize} + \log(2)) -
 * 1) \f$
 */
template <typename T> class ExpM1Axis : public Axis<T> {

  void print(std::ostream &os) const {
    os << "ExpM1Axis(";
    Axis<T>::print(os);
    os << ")";
  }

public:
  ExpM1Axis(T _low, T _high, T _stepsize = 1);

  ExpM1Axis(T _low, T _high, size_t _nodes);

  T transform(T x) const final;
  T back_transform(T t) const final;

  T derive(T x) const final;
  T back_derive(T t) const final;
};

/**
 * @brief Linear axis to describe data which varify not to much in orders of
 * magnitudes. It's the fastest axis evaluation.
 */
template <typename T> class LinAxis : public Axis<T> {
  void print(std::ostream &os) const {
    os << "LinAxis(";
    Axis<T>::print(os);
    os << ")";
  }

public:
  /**
   * @brief Linear Axis initalized with stepsize.
   */
  LinAxis(T _low, T _high, T _stepsize);

  /**
   * @brief Linear Axis with number of nodes.
   */
  LinAxis(T _low, T _high, size_t _nodes);

  T transform(T x) const final;
  T back_transform(T x) const final;
  T derive(T) const final;
  T back_derive(T x) const final;
};
} // namespace cubic_splines
