cmake_minimum_required(VERSION 3.14)

project(CubicInterpolation LANGUAGES CXX)

set(CubicInterpolation_VERSION_MAJOR 0)
set(CubicInterpolation_VERSION_MINOR 1)
set(CubicInterpolation_VERSION_PATCH 0)
set(CubicInterpolation_VERSION ${CubicInterpolation_VERSION_MAJOR}.${CubicInterpolation_VERSION_MINOR}.${CubicInterpolation_VERSION_PATCH})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

set(CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR} ${CMAKE_MODULE_PATH})
set(CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR} ${CMAKE_PREFIX_PATH})

add_subdirectory(src)
add_subdirectory(example)
# add_subdirectory(docs)

# if(BUILD_TESTING)
#     include(googletest)
#     enable_testing()
#     add_subdirectory(tests)
# endif()
