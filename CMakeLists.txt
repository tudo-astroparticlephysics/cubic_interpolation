cmake_minimum_required(VERSION 3.18.4)

project(test LANGUAGES CXX)

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

set(Boost_USE_DEBUG_RUNTIME OFF)
set(Boost_USE_STATIC_RUNTIME ON)
set(Boost_USE_DEBUG_LIBS OFF)
set(Boost_USE_STATIC_LIBS ON)
# Project wide setup

set(BOOST_LIBS serialization)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(eigen)
include(boost)

add_subdirectory(interpolate)

add_executable(test example.cxx)
target_link_libraries(test INTERPOLATE)
