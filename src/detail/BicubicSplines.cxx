#include "CubicInterpolation/BicubicSplines.h"
#include "CubicInterpolation/InterpolantBuilder.h"

#include <Eigen/Dense>
#include <boost/math/differentiation/finite_difference.hpp>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#include <boost/serialization/access.hpp>
#include <cmath>
#include <vector>

#include <chrono>
#include <iostream>

namespace cubic_splines {
namespace detail {

template <typename T> auto calculate_node(T val, unsigned int n_max) {
  auto n = std::floor(val);
  if (n < 0)
    return 0u;
  if (static_cast<unsigned int>(n) > n_max - 2u)
    return n_max - 2u;
  return static_cast<unsigned int>(n);
}

template <typename T> auto exponent_vector(T x) {
  auto x_vec = ::Eigen::Matrix<T, 4, 1>(4);
  for (size_t i = 0; i < 4; ++i) {
    x_vec(i) = 1;
    for (size_t j = 1; j <= i; ++j)
      x_vec(i) *= x;
  }
  return x_vec;
}
} // namespace detail

namespace detail {}

template <typename T> struct BicubicSplines<T>::RuntimeData {
  using MatrixX = ::Eigen::Matrix<T, ::Eigen::Dynamic, ::Eigen::Dynamic>;
  using Matrix4 = ::Eigen::Matrix<T, 4, 4>;

  MatrixX y, dydx1, dydx2, d2ydx1dx2;
  Matrix4 m =
      (Matrix4() << 1, 0, -3, 2, 0, 0, 3, -2, 0, 1, -2, 1, 0, 0, -1, 1).finished();

  RuntimeData() = default;

  RuntimeData(size_t n1, size_t n2)
      : y(n1, n2), dydx1(n1, n2), dydx2(n1, n2), d2ydx1dx2(n1, n2){};

  template <typename T1>
  RuntimeData(T1 &&_y, T1 &&_dydx1, T1 &&_dydx2, T1 &&_d2ydx1dx2)
      : y(std::forward<T1>(_y)), dydx1(std::forward<T1>(_dydx1)),
        dydx2(std::forward<T1>(_dydx2)), d2ydx1dx2(std::forward<T1>(_d2ydx1dx2)){};

  template <typename T1> static ::std::vector<T> to_vector(T1 m) {
    return ::std::vector<T>(m.data(), m.data() + m.rows() * m.cols());
  }

  auto get_dimensions() const {
    return std::array<long int, 2>{static_cast<long int>(y.rows()),
                                   static_cast<long int>(y.cols())};
  }

  StorageData to_storage_data() const;
};

template <typename T> struct BicubicSplines<T>::StorageData {
  using MatrixX = ::Eigen::Matrix<T, ::Eigen::Dynamic, ::Eigen::Dynamic>;

  ::std::array<long int, 2> size;
  ::std::vector<T> y, dydx1, dydx2, d2ydx1dx2;

  friend class boost::serialization::access;
  template <class Archive> void serialize(Archive &ar, const unsigned int) {
    ar &size;
    ar &y;
    ar &dydx1;
    ar &dydx2;
    ar &d2ydx1dx2;
  }

  template <typename V> inline auto to_matrix(V v) const {
    return ::Eigen::Map<MatrixX>(v.data(), size[0], size[1]);
  }

public:
  StorageData() = default;

  template <typename T1>
  StorageData(std::array<long int, 2> _size, T1 &&_y, T1 &&_dydx1, T1 &&_dydx2,
              T1 &&_d2ydx1dx2)
      : size(std::move(_size)), y(std::forward<T1>(_y)), dydx1(std::forward<T1>(_dydx1)),
        dydx2(std::forward<T1>(_dydx2)), d2ydx1dx2(std::forward<T1>(_d2ydx1dx2)){};

  auto to_runtime_data() const {
    return RuntimeData(to_matrix(y), to_matrix(dydx1), to_matrix(dydx2),
                       to_matrix(d2ydx1dx2));
  }
};

template <typename T>
typename BicubicSplines<T>::StorageData
BicubicSplines<T>::RuntimeData::to_storage_data() const {
  return StorageData(get_dimensions(), to_vector(y), to_vector(dydx1), to_vector(dydx2),
                     to_vector(d2ydx1dx2));
}

template <typename T>
BicubicSplines<T>::BicubicSplines(BicubicSplines::RuntimeData _data)
    : data(::std::make_unique<BicubicSplines::RuntimeData>(_data)) {}

template <typename T>
BicubicSplines<T>::BicubicSplines(Definition const &def, std::string path,
                                  std::string filename) {
  try {
    auto storage_data = load<BicubicSplines>(path, filename);
    *this = BicubicSplines(storage_data.to_runtime_data());
  } catch (std::system_error const &ex) {
    if (ex.code().value() != ENOENT)
      throw(ex);
    *this = BicubicSplines(def);
    save(data->to_storage_data(), path, filename);
  }
}

template <typename T>
std::tuple<T, T> BicubicSplines<T>::back_transform(Definition const &def,
                                                   unsigned long n1,
                                                   unsigned long n2) const {
  auto x1 = def.axis[0]->back_transform(n1);
  auto x2 = def.axis[1]->back_transform(n2);
  return std::make_tuple(x1, x2);
}

template <typename T>
template <typename T1>
std::array<T, 2>
BicubicSplines<T>::_prime(T1 func, std::array<std::unique_ptr<Axis<T>>, 2> const &axis,
                          unsigned int n1, unsigned int n2) {
  using boost::math::differentiation::finite_difference_derivative;
  auto f_x1 = [&func, ax = axis[0].get(), x2 = axis[1]->back_transform(n2)](T t1) {
    return func(ax->back_transform(t1), x2);
  };
  auto f_x2 = [&func, ax = axis[1].get(), x1 = axis[0]->back_transform(n1)](T t2) {
    return func(x1, ax->back_transform(t2));
  };
  auto dydx = std::array<T, 2>();
  dydx[0] = finite_difference_derivative(f_x1, static_cast<T>(n1));
  dydx[1] = finite_difference_derivative(f_x2, static_cast<T>(n2));
  return dydx;
}

/* template <typename T> */
/* template <typename T1> */
/* std::vector<std::tuple<unsigned int, T>> */
/* BicubicSplines<T>::_prime(std::vector<T> const &yi, T1 func, unsigned int n_max) { */
/*   auto spline = boost::math::interpolators::cardinal_cubic_b_spline<T>( */
/*       yi.data(), yi.size(), 0, 1, func(0u), func(n_max)); */
/*   auto diff = std::vector<std::tuple<unsigned int, T>>(); */
/*   for (auto i = 0u; i < n_max; ++i) */
/*     diff.emplace_back(i, spline.prime(i)); */
/*   return diff; */
/* } */

template <typename T>
template <typename T1>
T BicubicSplines<T>::_double_prime(Definition const &def, T1 func, unsigned int n1,
                                   unsigned int n2) {
  using boost::math::differentiation::finite_difference_derivative;
  auto f_x1 = [&func, ax = def.axis[0].get()](T t1, T t2) {
    return func(ax->back_transform(t1), t2);
  };
  auto dydx1 = [&f_x1, n1](T t2) {
    return finite_difference_derivative([&f_x1, t2](T t1) { return f_x1(t1, t2); },
                                        static_cast<T>(n1));
  };
  auto dydx1_x2 = [&dydx1, ax = def.axis[1].get()](T t2) {
    return dydx1(ax->back_transform(t2));
  };
  return finite_difference_derivative(dydx1_x2, static_cast<T>(n2));
}

template <typename T>
BicubicSplines<T>::BicubicSplines(Definition const &def)
    : data(std::make_shared<RuntimeData>(def.axis[0]->required_nodes(),
                                         def.axis[1]->required_nodes())) {
  using boost::math::differentiation::finite_difference_derivative;
  auto func = [&def](T x1, T x2) {
    if (def.f_trafo)
      return def.f_trafo->transform(def.f(x1, x2));
    return def.f(x1, x2);
  };

  for (auto n1 = 0u; n1 < data->y.rows(); ++n1) {
    for (auto n2 = 0u; n2 < data->y.cols(); ++n2) {
      auto x = back_transform(def, n1, n2);
      data->y(n1, n2) = func(std::get<0>(x), std::get<1>(x));
    }
  }

  auto n_rows = data->y.rows();
  auto n_cols = data->y.cols();
  if (def.approx_derivates) {
    for (auto n1 = 0u; n1 < n_cols; ++n1) {
      auto yi = RuntimeData::to_vector(data->y.col(n1));
      auto diff = [this, &def, &func, n1](unsigned int n) {
        return _prime(func, def.axis, n, n1)[0];
      };
      auto spline = boost::math::interpolators::cardinal_cubic_b_spline<T>(
          yi.data(), yi.size(), 0, 1, diff(0u), diff(n_rows - 1));
      for (auto row = 0; row < n_rows; ++row)
        data->dydx1(n1, row) = spline.prime(row);
    }
    data->dydx1.transposeInPlace();

    ::Eigen::Matrix<T, ::Eigen::Dynamic, ::Eigen::Dynamic> y_rowise = data->y;
    y_rowise.transposeInPlace();
    for (auto n2 = 0u; n2 < n_rows; ++n2) {
      auto yi = RuntimeData::to_vector(y_rowise.col(n2));
      auto diff = [this, &def, &func, n2](unsigned int n) {
        return _prime(func, def.axis, n2, n)[1];
      };
      auto spline = boost::math::interpolators::cardinal_cubic_b_spline<T>(
          yi.data(), yi.size(), 0, 1, diff(0u), diff(n_cols - 1));
      for (auto col = 0; col < n_cols; ++col)
        data->dydx2(col, n2) = spline.prime(col);
    }
    data->dydx2.transposeInPlace();
  } else {
    for (auto n1 = 0u; n1 < n_rows; ++n1) {
      for (auto n2 = 0u; n2 < n_cols; ++n2) {
        auto dydx = _prime(func, def.axis, n1, n2);
        data->dydx1(n1, n2) = dydx[0];
        data->dydx2(n1, n2) = dydx[1];
      }
    }
  }
  if (def.approx_derivates) {
    ::Eigen::Matrix<T, ::Eigen::Dynamic, ::Eigen::Dynamic> dydx1_rowise = data->dydx1;
    dydx1_rowise.transposeInPlace();
    for (auto n2 = 0u; n2 < n_rows; ++n2) {
      auto yi = RuntimeData::to_vector(dydx1_rowise.col(n2));
      auto diff = [this, &def, &func, n2](unsigned int n) {
        return _double_prime(def, func, n2, n);
      };
      auto spline = boost::math::interpolators::cardinal_cubic_b_spline<T>(
          yi.data(), yi.size(), 0, 1, diff(0u), diff(n_cols - 1));
      for (auto col = 0; col < n_cols; ++col)
        data->d2ydx1dx2(col, n2) = spline.prime(col);
    }
    data->d2ydx1dx2.transposeInPlace();
  } else {
    for (auto n1 = 0u; n1 < n_rows; ++n1) {
      for (auto n2 = 0u; n2 < n_cols; ++n2)
        data->d2ydx1dx2(n1, n2) = _double_prime(def, func, n1, n2);
    }
  }
}

template <typename T> T BicubicSplines<T>::evaluate(T x0, T x1) const {
  auto n0 = detail::calculate_node(x0, data->y.cols());
  auto n1 = detail::calculate_node(x1, data->y.rows());
  auto temp = ::Eigen::Matrix<T, 4, 4>(4, 4);
  temp.template block<2, 2>(0, 0) = data->y.block(n0, n1, 2, 2);
  temp.template block<2, 2>(2, 0) = data->dydx1.block(n0, n1, 2, 2);
  temp.template block<2, 2>(0, 2) = data->dydx2.block(n0, n1, 2, 2);
  temp.template block<2, 2>(2, 2) = data->d2ydx1dx2.block(n0, n1, 2, 2);
  auto v1 = detail::exponent_vector(x0 - n0);
  auto v2 = detail::exponent_vector(x1 - n1);
  return v1.dot((data->m.transpose() * (temp * data->m)) * v2);
}

template <typename T> std::array<T, 2> BicubicSplines<T>::prime(T x0, T x1) const {
  using boost::math::differentiation::finite_difference_derivative;
  auto grad = std::array<T, 2>();
  grad[0] =
      finite_difference_derivative([this, x1](T x_0) { return evaluate(x_0, x1); }, x0);
  grad[1] =
      finite_difference_derivative([this, x0](T x_1) { return evaluate(x0, x_1); }, x1);
  return grad;
}

template <typename T> T BicubicSplines<T>::double_prime(T x0, T x1) const {
  using boost::math::differentiation::finite_difference_derivative;
  return finite_difference_derivative(
      [this, x1](T x_1) {
        return finite_difference_derivative(
            [this, x_1](T x_2) { return evaluate(x_1, x_2); }, x1);
      },
      x0);
}
} // namespace cubic_splines

template class cubic_splines::BicubicSplines<double>;
template class cubic_splines::BicubicSplines<float>;
