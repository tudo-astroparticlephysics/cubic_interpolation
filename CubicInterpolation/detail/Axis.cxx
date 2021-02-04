#include "../Axis.h"
#include <cmath>

using namespace cubic_splines;

template <typename T>
Axis<T>::Axis(T _low, T _high, T _stepsize)
    : low(_low), high(_high), stepsize(_stepsize) {}

template <typename T> unsigned int Axis<T>::required_nodes() const {
  auto nodes = transform(high);
  if (std::floor(nodes) == nodes)
    return nodes + 1;
  return nodes + 2;
}

template <typename T>
ExpAxis<T>::ExpAxis(T _low, T _high, T _stepsize) : Axis<T>(_low, _high, _stepsize) {}

template <typename T>
ExpAxis<T>::ExpAxis(T _low, T _high, size_t _nodes) : Axis<T>(_low, _high, 0.f) {
  this->stepsize = std::log(_high / _low) / (_nodes - 1);
}

template <typename T> T ExpAxis<T>::transform(T x) const {
  return std::log(x / this->low) / this->stepsize;
}

template <typename T> T ExpAxis<T>::back_transform(T t) const {
  return this->low * std::exp(t * this->stepsize);
}
template <typename T> T ExpAxis<T>::derive(T x) const {
  return 1. / (x * this->stepsize);
}

template <typename T> T ExpAxis<T>::back_derive(T t) const {
  return this->low * this->stepsize * std::exp(t * this->stepsize);
}

template <typename T>
LinAxis<T>::LinAxis(T _low, T _high, T _stepsize) : Axis<T>(_low, _high, _stepsize) {}

template <typename T>
LinAxis<T>::LinAxis(T _low, T _high, size_t _nodes) : Axis<T>(_low, _high, 0.f) {
  this->stepsize = (_high - _low) / (_nodes - 1);
}

template <typename T> T LinAxis<T>::transform(T x) const {
  return (x - this->low) / this->stepsize;
}

template <typename T> T LinAxis<T>::back_transform(T x) const {
  return x * this->stepsize + this->low;
}

template <typename T> T LinAxis<T>::derive(T) const { return 1. / this->stepsize; }
template <typename T> T LinAxis<T>::back_derive(T) const { return this->stepsize; }

namespace cubic_splines {
template class Axis<double>;
template class Axis<float>;
template class LinAxis<double>;
template class LinAxis<float>;
template class ExpAxis<double>;
template class ExpAxis<float>;
} // namespace cubic_splines
