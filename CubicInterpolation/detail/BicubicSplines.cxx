#include "../BicubicSplines.h"
#include "../InterpolantBuilder.h"

#include <Eigen/Dense>
#include <boost/math/differentiation/finite_difference.hpp>
#include <boost/serialization/access.hpp>
#include <cmath>
#include <vector>

namespace cubic_splines {
namespace detail {
unsigned int calculate_node(float val, unsigned int n_max) {
  auto n = std::floor(val);
  if (n < 0) {
    return 0;
  } else if (n > n_max - 2u) {
    return n_max - 2u;
  }
  return n;
}

auto exponent_vector(float x) {
  auto x_vec = ::Eigen::Vector4f(4);
  for (size_t i = 0; i < 4; ++i) {
    x_vec(i) = 1;
    for (size_t j = 1; j <= i; ++j)
      x_vec(i) *= x;
  }
  return x_vec;
}
} // namespace detail

struct BicubicSplines::RuntimeData {
  ::Eigen::MatrixXf y, dydx1, dydx2, d2ydx1dx2;
  ::Eigen::Matrix4f m = (::Eigen::Matrix4f() << 1, 0, -3, 2, 0, 0, 3, -2, 0, 1,
                         -2, 1, 0, 0, -1, 1)
                            .finished();

  RuntimeData() = default;

  RuntimeData(size_t n1, size_t n2)
      : y(n1, n2), dydx1(n1, n2), dydx2(n1, n2), d2ydx1dx2(n1, n2){};

  RuntimeData(::Eigen::MatrixXf _y, ::Eigen::MatrixXf _dydx1,
              ::Eigen::MatrixXf _dydx2, ::Eigen::MatrixXf _d2ydx1dx2)
      : y(std::move(_y)), dydx1(std::move(_dydx1)), dydx2(std::move(_dydx2)),
        d2ydx1dx2(std::move(_d2ydx1dx2)){};

  StorageData to_storage_data() const;
};

class BicubicSplines::StorageData {
  ::std::array<size_t, 2> size;
  ::std::vector<float> y, dydx1, dydx2, d2ydx1dx2;

  friend class boost::serialization::access;
  template <class Archive> void serialize(Archive &ar, const unsigned int) {
    ar &size;
    ar &y;
    ar &dydx1;
    ar &dydx2;
    ar &d2ydx1dx2;
  }

  template <typename M, typename S> inline auto to_vector(M m, S s) {
    return ::std::vector<float>(m.data(), m.data() + s[0] * s[1]);
  }

  template <typename V, typename S> inline auto to_matrix(V v, S s) {
    return ::Eigen::Map<::Eigen::MatrixXf>(v.data(), s[0], s[1]);
  }

public:
  StorageData() = default;
  using matrix_t =
      ::std::decay<decltype(BicubicSplines::RuntimeData().y)>::type;
  StorageData(matrix_t _y, matrix_t _dydx1, matrix_t _dydx2,
              matrix_t _d2ydx1dx2)
      : size({(size_t)_y.rows(), (size_t)_y.cols()}), y(to_vector(_y, size)),
        dydx1(to_vector(_dydx1, size)), dydx2(to_vector(_dydx2, size)),
        d2ydx1dx2(to_vector(_d2ydx1dx2, size)){};

  auto to_runtime_data() {
    return RuntimeData{to_matrix(y, size), to_matrix(dydx1, size),
                       to_matrix(dydx2, size), to_matrix(d2ydx1dx2, size)};
  }
};

BicubicSplines::StorageData
BicubicSplines::RuntimeData::to_storage_data() const {
  return StorageData(y, dydx1, dydx2, d2ydx1dx2);
}

/* BicubicSplines &BicubicSplines::operator=(BicubicSplines const &splines) { */
/*   data = std::make_unique<RuntimeData>(*splines.data); */
/*   return *this; */
/* } */

BicubicSplines::BicubicSplines(BicubicSplines::RuntimeData _data)
    : data(::std::make_unique<BicubicSplines::RuntimeData>(_data)) {}

BicubicSplines::BicubicSplines(Definition const &def, std::string path,
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

BicubicSplines::BicubicSplines(Definition const &def)
    : data(std::make_unique<RuntimeData>(def.axis[0]->required_nodes(),
                                         def.axis[1]->required_nodes())) {
  using boost::math::differentiation::finite_difference_derivative;
  auto func = [&def](float x1, float x2) {
    if (def.f_trafo)
      return def.f_trafo->transform(def.f(x1, x2));
    return def.f(x1, x2);
  };
  for (auto n1 = 0u; n1 < data->y.rows(); ++n1) {
    for (auto n2 = 0u; n2 < data->y.cols(); ++n2) {
      auto x1 = def.axis[0]->back_transform(n1);
      auto x2 = def.axis[1]->back_transform(n2);
      auto dfdx1 = def.axis[0]->derive(x1);
      auto dfdx2 = def.axis[1]->derive(x2);
      data->y(n1, n2) = func(x1, x2);
      data->dydx1(n1, n2) =
          finite_difference_derivative(
              [this, &func, x2](float x) { return func(x, x2); }, x1) *
          dfdx1;
      data->dydx2(n1, n2) =
          finite_difference_derivative(
              [this, &func, x1](float x) { return func(x1, x); }, x2) *
          dfdx2;
      data->d2ydx1dx2(n1, n2) =
          finite_difference_derivative(
              [this, &func, x1, x2, dfdx1, dfdx2](float x_1) {
                return finite_difference_derivative(
                           [this, &func, x_1, dfdx2](float x_2) {
                             return func(x_1, x_2);
                           },
                           x2) *
                       dfdx2;
              },
              x1) *
          dfdx1;
    }
  }
}

/* BicubicSplines::~BicubicSplines() {} */

float BicubicSplines::evaluate(float x0, float x1) const {
  auto n0 = detail::calculate_node(x0, data->y.cols());
  auto n1 = detail::calculate_node(x1, data->y.rows());
  auto temp = Eigen::Matrix4f(4, 4);
  temp.block<2, 2>(0, 0) = data->y.block(n0, n1, 2, 2);
  temp.block<2, 2>(2, 0) = data->dydx1.block(n0, n1, 2, 2);
  temp.block<2, 2>(0, 2) = data->dydx2.block(n0, n1, 2, 2);
  temp.block<2, 2>(2, 2) = data->d2ydx1dx2.block(n0, n1, 2, 2);
  auto v1 = detail::exponent_vector(x0 - n0);
  auto v2 = detail::exponent_vector(x1 - n1);
  return v1.dot((data->m.transpose() * (temp * data->m)) * v2);
}

std::array<float, 2> BicubicSplines::prime(float x0, float x1) const {
  using boost::math::differentiation::finite_difference_derivative;
  auto grad = std::array<float, 2>();
  grad[0] = finite_difference_derivative(
      [this, x1](float x_0) { return evaluate(x_0, x1); }, x0);
  grad[1] = finite_difference_derivative(
      [this, x0](float x_1) { return evaluate(x0, x_1); }, x1);
  return grad;
}

float BicubicSplines::double_prime(float x0, float x1) const {
  using boost::math::differentiation::finite_difference_derivative;
  return finite_difference_derivative(
      [this, x0, x1](float x_1) {
        return finite_difference_derivative(
            [this, x1, x_1](float x_2) { return evaluate(x_1, x_2); }, x1);
      },
      x0);
}
} // namespace cubic_splines
