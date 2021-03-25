cmake_minimum_required(VERSION 3.12)

project(CubicInterpolation LANGUAGES CXX)

if (EXISTS ${CMAKE_BINARY_DIR}/conan_paths.cmake)
    include(${CMAKE_BINARY_DIR}/conan_paths.cmake)
endif()

message(STATUS ${CMAKE_PREFIX_PATH})
message(STATUS ${CMAKE_MODULE_PATH})

set(CubicInterpolation_VERSION_MAJOR 0)
set(CubicInterpolation_VERSION_MINOR 1)
set(CubicInterpolation_VERSION_PATCH 3)
set(CubicInterpolation_VERSION ${CubicInterpolation_VERSION_MAJOR}.${CubicInterpolation_VERSION_MINOR}.${CubicInterpolation_VERSION_PATCH})

set(CMAKE_CXX_STANDARD 14)

find_package(Eigen3 REQUIRED)
find_package(Boost COMPONENTS filesystem serialization REQUIRED)

add_subdirectory(src)

option(BUILD_EXAMPLE "build example" OFF)
if(BUILD_EXAMPLE)
    add_subdirectory(example)
endif()

option(BUILD_DOCUMENTATION "build documentation" OFF)
if(BUILD_DOCUMENTATION)
    add_subdirectory(docs)
endif()

if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()
