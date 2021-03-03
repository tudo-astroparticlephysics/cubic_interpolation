cmake_minimum_required(VERSION 3.12)

project(CubicInterpolation LANGUAGES CXX)

set(CubicInterpolation_VERSION_MAJOR 0)
set(CubicInterpolation_VERSION_MINOR 1)
set(CubicInterpolation_VERSION_PATCH 3)
set(CubicInterpolation_VERSION ${CubicInterpolation_VERSION_MAJOR}.${CubicInterpolation_VERSION_MINOR}.${CubicInterpolation_VERSION_PATCH})

set(CMAKE_CXX_STANDARD 14)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup(TARGETS)

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
