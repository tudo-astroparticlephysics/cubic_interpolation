cmake_minimum_required(VERSION 3.14)

project(test LANGUAGES CXX)

include(GNUInstallDirs)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

add_subdirectory(vendor)
add_subdirectory(interpolate)

add_executable(test example.cxx)
target_link_libraries(test INTERPOLATE)
