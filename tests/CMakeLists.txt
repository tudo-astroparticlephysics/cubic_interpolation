enable_testing()

include(GoogleTest)

find_package(GTest REQUIRED)

set(Libs
    CubicInterpolation::CubicInterpolation
    GTest::GTest)

# add_executable(TestAxis TestAxis.cpp)
# target_link_libraries(TestAxis PRIVATE gtest CubicInterpolation)
# gtest_discover_tests(TestAxis)

add_executable(TestCubicSplines TestCubicSplines.cpp)
target_link_libraries(TestCubicSplines ${Libs})
gtest_discover_tests(TestCubicSplines)

add_executable(TestBicubicSplines TestBicubicSplines.cpp)
target_link_libraries(TestBicubicSplines ${Libs})
gtest_discover_tests(TestBicubicSplines)

add_executable(TestFindParameter TestFindParameter.cpp)
target_link_libraries(TestFindParameter ${Libs})
gtest_discover_tests(TestFindParameter)
