cmake_minimum_required(VERSION 3.14)

project(test LANGUAGES CXX)

include(GNUInstallDirs)
include(CTest)

set(project_BUILD_TESTING BUILD_TESTING)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# set(Boost_USE_DEBUG_RUNTIME FALSE)
set(Boost_USE_STATIC_RUNTIME OFF)
# set(Boost_USE_DEBUG_LIBS OFF)
set(Boost_USE_STLPORT OFF)
set(Boost_USE_STATIC_LIBS OFF)

message(STATUS "BUILD_SHARED_LIBS: ${BUILD_SHARED_LIBS}")
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(eigen)

set(BOOST_LIBS serialization filesystem)
include(boost)

add_subdirectory(interpolate)

add_subdirectory(example)

if(project_BUILD_TESTING)
    include(googletest)

    enable_testing()
    add_subdirectory(tests)
endif()
