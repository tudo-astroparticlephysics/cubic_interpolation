include(GoogleTest)

add_executable(TestAxis TestAxis.cpp)
target_link_libraries(TestAxis PRIVATE gtest)
gtest_discover_tests(TestAxis)
