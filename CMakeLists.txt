cmake_minimum_required(VERSION 3.13)

project(CubicInterpolation LANGUAGES CXX)

set(CubicInterpolation_VERSION_MAJOR 0)
set(CubicInterpolation_VERSION_MINOR 1)
set(CubicInterpolation_VERSION_PATCH 0)
set(CubicInterpolation_VERSION ${CubicInterpolation_VERSION_MAJOR}.${CubicInterpolation_VERSION_MINOR}.${CubicInterpolation_VERSION_PATCH})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

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
