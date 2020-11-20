include(GoogleTest)

add_executable(TestAxis TestAxis.cpp)
target_link_libraries(TestAxis PRIVATE gtest CubicInterpolation)
gtest_discover_tests(TestAxis)

add_executable(TestCubicSplines TestCubicSplines.cpp)
target_link_libraries(TestCubicSplines PRIVATE gtest CubicInterpolation)
gtest_discover_tests(TestCubicSplines)
