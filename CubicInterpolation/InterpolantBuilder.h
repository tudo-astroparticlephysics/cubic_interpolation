#pragma once

#include <cerrno>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <system_error>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>

#include <boost/filesystem.hpp>

namespace cubic_splines {

namespace fs = ::boost::filesystem;

template <typename T1, typename T2 = typename T1::StorageData>
T2 load(fs::path path, fs::path filename) {
  if (fs::is_regular_file(path / filename)) {
    std::ifstream ifs((path / filename).c_str());
    auto Data = T2();
    if (ifs.is_open()) {
      boost::archive::text_iarchive ia(ifs);
      ia >> Data;
      return Data;
    }
  }
  throw std::system_error(ENOENT, std::generic_category(),
                          "Interpolation tables couldn't be found");
};

template <typename T>
bool save(T const &storage_data, fs::path path, fs::path filename) {
  if (not fs::exists(path / filename)) {
    std::ofstream ofs((path / filename).c_str());
    while (ofs.good()) {
      boost::archive::text_oarchive oa(ofs);
      oa << storage_data;
      return true;
    }
  }
  return false;
};

/* template <typename T1, typename T2 = typename T1::Definition, */
/*           typename T3 = typename T1::StorageData> */
/* class InterpolantBuilder { */

/*   T1 load(boost::filesystem::path path, boost::filesystem::path filename) {
 */
/*     if (boost::filesystem::is_regular_file(path / filename)) { */
/*       std::ifstream ifs((path / filename).c_str()); */
/*       auto Data = T3(); */
/*       if (ifs.is_open()) { */
/*         boost::archive::text_iarchive ia(ifs); */
/*         ia >> Data; */
/*         return T1(Data); */
/*       } */
/*     } */
/*     throw std::system_error(ENOENT, std::generic_category(), */
/*                             "Interpolation tables couldn't be found"); */
/*   }; */

/*   template <typename... T> */
/*   bool save(boost::filesystem::path path, boost::filesystem::path filename,
 */
/*             T... args) { */
/*     if (not boost::filesystem::exists(path / filename)) { */
/*       auto data = T3(args...); */
/*       std::ofstream ofs((path / filename).c_str()); */
/*       while (ofs.good()) { */
/*         boost::archive::text_oarchive oa(ofs); */
/*         oa << data; */
/*         return true; */
/*       } */
/*     } */
/*     return false; */
/*   }; */

/*   template <typename... T> auto transform(T2 const &def, T... args) { */
/*     auto fx = def.f(args...); */
/*     if (def.f_trafo) */
/*       fx = def.f_trafo->transform(fx); */
/*     return fx; */
/*   }; */

/* public: */
/*   InterpolantBuilder() = default; */

/* T1 build(T2 const &def, std::string save_path, std::string filename); */
/* }; */

/* template <> */
/* BicubicSplines */
/* InterpolantBuilder<BicubicSplines>::build(BicubicSplines::Definition const &,
 */
/*                                           std::string, std::string); */

/* template <> */
/* CubicSplines */
/* InterpolantBuilder<CubicSplines>::build(CubicSplines::Definition const &, */
/*                                         std::string, std::string); */

} // namespace cubic_splines
