cmake_minimum_required(VERSION 3.1)
project(CubicInterpolationTest CXX)

set(CMAKE_CXX_STANDARD 17)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

enable_testing()

include(GoogleTest)

# add_executable(TestAxis TestAxis.cpp)
# target_link_libraries(TestAxis PRIVATE gtest CubicInterpolation)
# gtest_discover_tests(TestAxis)

add_executable(TestCubicSplines TestCubicSplines.cpp)
target_link_libraries(TestCubicSplines ${CONAN_LIBS})
gtest_discover_tests(TestCubicSplines)

add_executable(TestBicubicSplines TestBicubicSplines.cpp)
target_link_libraries(TestBicubicSplines ${CONAN_LIBS})
gtest_discover_tests(TestBicubicSplines)

add_executable(TestFindParameter TestFindParameter.cpp)
target_link_libraries(TestFindParameter ${CONAN_LIBS})
gtest_discover_tests(TestFindParameter)
