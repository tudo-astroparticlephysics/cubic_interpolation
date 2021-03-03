#pragma once

#include <cerrno>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <system_error>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
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
      boost::archive::binary_iarchive ia(ifs);
      ia >> Data;
      return Data;
    }
  }
  throw std::system_error(ENOENT, std::generic_category(),
                          "Interpolation tables couldn't be found");
};

template <typename T> bool save(T const &storage_data, fs::path path, fs::path filename) {
  if (!fs::exists(path / filename)) {
    std::ofstream ofs((path / filename).c_str());
    while (ofs.good()) {
      boost::archive::binary_oarchive oa(ofs);
      oa << storage_data;
      return true;
    }
  }
  return false;
};
} // namespace cubic_splines
