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
struct CubicSplines::StorageData {
  std::vector<float> y;
  float lower_lim_derivate;
  float upper_lim_derivate;

  friend class boost::serialization::access;
  template <class Archive> void serialize(Archive &ar, const unsigned int) { ar &y; };

public:
  StorageData() = default;

  template <typename T>
  StorageData(T const &_y, float _lower_lim_derivate, float _upper_lim_derivate)
      : y(_y.begin(), _y.end()), lower_lim_derivate(_lower_lim_derivate),
        upper_lim_derivate(_upper_lim_derivate){};

  RuntimeData to_runtime_data() const;
};

struct CubicSplines::RuntimeData : public CubicSplines::StorageData {

  boost::math::interpolators::cardinal_cubic_b_spline<float> spline;

  RuntimeData() = default;

  template <typename T>
  RuntimeData(T _data, float _lower_lim_derivate, float _upper_lim_derivate)
      : StorageData(_data, _lower_lim_derivate, _upper_lim_derivate),
        spline(y.data(), y.size(), 0, 1, lower_lim_derivate, upper_lim_derivate) {}

  StorageData to_storage_data() const { return StorageData(*this); };
};

CubicSplines::RuntimeData CubicSplines::StorageData::to_runtime_data() const {
  return RuntimeData(y, lower_lim_derivate, upper_lim_derivate);
};

CubicSplines::CubicSplines(CubicSplines::RuntimeData _data)
    : data(::std::make_unique<CubicSplines::RuntimeData>(_data)) {}

CubicSplines::CubicSplines(Definition const &def, std::string path,
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

CubicSplines::CubicSplines(Definition const &def)
    : data(::std::make_unique<CubicSplines::RuntimeData>()) {
  using boost::math::differentiation::finite_difference_derivative;
  auto y = std::vector<float>(def.axis->required_nodes());
  auto func = [&def](float x) {
    auto fx = def.f(x);
    if (def.f_trafo)
      fx = def.f_trafo->transform(fx);
    return fx;
  };
  for (size_t n = 0; n < y.size(); ++n)
    y[n] = func(def.axis->back_transform(n));
  auto low = def.axis->back_transform(0);
  auto up = def.axis->back_transform(y.size() - 1);
  auto lower_lim_derivate =
      finite_difference_derivative(func, low) * def.axis->derive(low);
  auto upper_lim_derivate =
      finite_difference_derivative(func, def.axis->back_transform(y.size() - 1)) *
      def.axis->derive(up);
  data = ::std::make_unique<CubicSplines::RuntimeData>(y, lower_lim_derivate,
                                                       upper_lim_derivate);
}

float CubicSplines::evaluate(float x) const { return data->spline(x); };

float CubicSplines::prime(float x) const { return data->spline.prime(x); };

float CubicSplines::double_prime(float x) const { return data->spline.double_prime(x); };

} // namespace cubic_splines
