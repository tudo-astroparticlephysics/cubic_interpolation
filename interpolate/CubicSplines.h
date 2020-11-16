#pragma once

#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#include <boost/serialization/access.hpp>
#include <cmath>
#include <functional>
#include <memory>
#include <vector>

#include "Axis.h"

class CubicSplines {
  boost::math::interpolators::cardinal_cubic_b_spline<float> spline;

public:
  CubicSplines() = default;

  template <typename T> CubicSplines(T val);

  static constexpr size_t N = 1;

  struct Definition;

  struct Data;

  float evaluate(int node, float x);
};

struct CubicSplines::Definition {
  std::function<float(float)> f;
  std::unique_ptr<Axis> f_trafo = nullptr;
  std::unique_ptr<Axis> axis;
};

class CubicSplines::Data {
  std::vector<float> y;

  friend class boost::serialization::access;
  template <class Archive> void serialize(Archive &ar, const unsigned int) {
    ar &y;
  }

public:
  Data() = default;

  template <typename T>
  Data(T&& _y) : y(std::forward<T>(_y)) {};

  CubicSplines build();
};

template <typename T>
CubicSplines::CubicSplines(T val) : spline(val.data(), val.size(), 0, 1) {}

float CubicSplines::evaluate(int node, float x) { return spline(x + node); }

CubicSplines CubicSplines::Data::build() {
    return CubicSplines(y);
}
