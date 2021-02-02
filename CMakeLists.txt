cmake_minimum_required(VERSION 3.14)

project(CubicInterpolation LANGUAGES CXX)

include(GNUInstallDirs)
include(CTest)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# set(Boost_USE_DEBUG_RUNTIME FALSE)
set(Boost_USE_STATIC_RUNTIME OFF)
# set(Boost_USE_DEBUG_LIBS OFF)
set(Boost_USE_STLPORT OFF)
set(Boost_USE_STATIC_LIBS OFF)
set(BOOST_LIBS serialization filesystem)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(eigen)
include(boost)

add_library(CubicInterpolation)
add_library(CubicInterpolation::CubicInterpolation ALIAS CubicInterpolation)

add_subdirectory(CubicInterpolation)

target_include_directories(CubicInterpolation INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/CubicInterpolation>
    )


set_target_properties(CubicInterpolation PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    )

install(TARGETS CubicInterpolation EXPORT CubicInterpolationTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )

install(EXPORT CubicInterpolationTargets
    FILE CubicInterpolationTargets.cmake
    NAMESPACE CubicInterpolation::
    DESTINATION lib/cmake/CubicInterpolation
    )

include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/CubicInterpolationConfig.cmake
    INSTALL_DESTINATION lib/cmake
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
    )

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/CubicInterpolationConfig.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/CubicInterpolation
    )

export(EXPORT CubicInterpolationTargets
    FILE "${CMAKE_CURRENT_BINARY_DIR}/CubicInterpolationTargets.cmake"
    )

add_subdirectory(example)
add_subdirectory(docs)

if(BUILD_TESTING)
    include(googletest)

    enable_testing()
    add_subdirectory(tests)
endif()
