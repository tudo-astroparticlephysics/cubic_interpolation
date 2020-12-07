#include "../Axis.h"

#include <cmath>

namespace cubic_splines {

void Axis::print(std::ostream &os) const {
  os << "low: " << low << ", high: " << high << ", stepsize: " << stepsize;
};

Axis::Axis(float _low, float _high, float _stepsize)
    : low(_low), high(_high), stepsize(_stepsize) {}

size_t Axis::required_nodes() const {
  auto nodes = transform(high);
  if (std::floor(nodes) == nodes)
    return nodes + 1;
  return nodes + 2;
}

std::ostream &operator<<(std::ostream &out, const Axis &axis) {
  axis.print(out);
  return out;
}

void ExpAxis::print(std::ostream &os) const {
  os << "ExpAxis(";
  Axis::print(os);
  os << ")";
}

ExpAxis::ExpAxis(float _low, float _high, float _stepsize)
    : Axis(_low, _high, _stepsize) {}

/**
 * @brief Exponential Axis initialized with number of nodes.
 */
ExpAxis::ExpAxis(float _low, float _high, size_t _nodes)
    : Axis(_low, _high, 0.f) {
  stepsize = std::log(high / low) / (_nodes - 1);
}

float ExpAxis::transform(float x) const noexcept {
  return std::log(x / low) / stepsize;
}

float ExpAxis::back_transform(float x) const noexcept {
  return low * std::exp(x * stepsize);
}

float ExpAxis::derive(float x) const { return (x * stepsize); }

void LinAxis::print(std::ostream &os) const {
  os << "LinAxis(";
  Axis::print(os);
  os << ")";
}

LinAxis::LinAxis(float _low, float _high, float _stepsize)
    : Axis(_low, _high, _stepsize) {}

LinAxis::LinAxis(float _low, float _high, size_t _nodes)
    : Axis(_low, _high, 0.f) {
  stepsize = (high - low) / (_nodes - 1);
}

float LinAxis::transform(float x) const noexcept {
  return (x - low) / stepsize;
}

float LinAxis::back_transform(float x) const noexcept {
  return x * stepsize + low;
}

float LinAxis::derive(float) const { return stepsize; }

} // namespace cubic_splines
