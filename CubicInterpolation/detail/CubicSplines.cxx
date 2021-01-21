#include "../CubicSplines.h"
#include "../InterpolantBuilder.h"

#include <boost/math/differentiation/finite_difference.hpp>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#include <boost/serialization/access.hpp>
#include <vector>

namespace cubic_splines {

/**
 * @brief Storage class to write and load the interpolation tables from  disk.
 * After reading and writing the object will be destructed.
 */
template <typename T> struct CubicSplines<T>::StorageData {
  std::vector<T> y;
  T lower_lim_derivate;
  T upper_lim_derivate;

  friend class boost::serialization::access;
  template <class Archive> void serialize(Archive &ar, const unsigned int) { ar &y; };

public:
  StorageData() = default;

  template <typename T1>
  StorageData(T1 const &_y, T _lower_lim_derivate, T _upper_lim_derivate)
      : y(_y.begin(), _y.end()), lower_lim_derivate(_lower_lim_derivate),
        upper_lim_derivate(_upper_lim_derivate){};

  auto to_runtime_data() const {
    return RuntimeData(y, lower_lim_derivate, upper_lim_derivate);
  }
};

template <typename T>
struct CubicSplines<T>::RuntimeData : public CubicSplines<T>::StorageData {

  boost::math::interpolators::cardinal_cubic_b_spline<T> spline;

  RuntimeData() = default;

  template <typename T1>
  RuntimeData(T1 _data, T _lower_lim_derivate, T _upper_lim_derivate)
      : StorageData(_data, _lower_lim_derivate, _upper_lim_derivate),
        spline(_data.data(), _data.size(), 0, 1, _lower_lim_derivate,
               _upper_lim_derivate) {}

  auto to_storage_data() const { return StorageData(*this); };
};

template <typename T>
CubicSplines<T>::CubicSplines(CubicSplines::RuntimeData _data)
    : data(::std::make_unique<CubicSplines::RuntimeData>(_data)) {}

template <typename T>
CubicSplines<T>::CubicSplines(Definition const &def, std::string path,
                              std::string filename) {
  try {
    auto storage_data = load<CubicSplines>(path, filename);
    *this = CubicSplines(storage_data.to_runtime_data());
  } catch (std::system_error const &ex) {
    if (ex.code().value() != ENOENT)
      throw(ex);
    *this = CubicSplines(def);
    save(data->to_storage_data(), path, filename);
  }
}

template <typename T>
CubicSplines<T>::CubicSplines(Definition const &def)
    : data(::std::make_unique<CubicSplines::RuntimeData>()) {
  using boost::math::differentiation::finite_difference_derivative;
  auto y = std::vector<T>(def.axis->required_nodes());
  auto func = [&def](T x) {
    auto fx = def.f(x);
    if (def.f_trafo)
      fx = def.f_trafo->transform(fx);
    return fx;
  };
  for (size_t n = 0; n < y.size(); ++n)
    y[n] = func(def.axis->back_transform(n));
  auto n_max = y.size() - 1;

  auto low = def.axis->back_transform(0);
  auto up = def.axis->back_transform(n_max);

  /* std::cout << "low: " << low << std::endl; */
  /* std::cout << "up: " << up << std::endl; */

  /* std::cout << "diff f: " << finite_difference_derivative(func, low) << std::endl; */
  /* std::cout << "diff f: " << finite_difference_derivative(func, up) << std::endl; */

  /* std::cout << "diff trafo: " << def.axis->back_derive(0) << std::endl; */
  /* std::cout << "diff trafo: " << def.axis->back_derive(n_max) << std::endl; */

  auto lower_lim_derivate =
      finite_difference_derivative(func, low) * def.axis->back_derive(0);
  auto upper_lim_derivate =
      finite_difference_derivative(func, up) * def.axis->back_derive(n_max);

  data = ::std::make_unique<CubicSplines::RuntimeData>(y, lower_lim_derivate,
                                                       upper_lim_derivate);
}

template <typename T> T CubicSplines<T>::evaluate(T x) const { return data->spline(x); };

template <typename T> T CubicSplines<T>::prime(T x) const {
  return data->spline.prime(x);
};

template <typename T> T CubicSplines<T>::double_prime(T x) const {
  return data->spline.double_prime(x);
};

} // namespace cubic_splines
