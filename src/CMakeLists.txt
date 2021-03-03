include(GNUInstallDirs)

add_library(CubicInterpolation)
add_library(CubicInterpolation::CubicInterpolation ALIAS CubicInterpolation)

if(MSVC)
    target_compile_definitions(CubicInterpolation PRIVATE _USE_MATH_DEFINES)
endif()

add_subdirectory(CubicInterpolation)
add_subdirectory(detail)

target_link_libraries(CubicInterpolation
    CONAN_PKG::eigen
    CONAN_PKG::boost
    )

target_include_directories(CubicInterpolation PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )

install(TARGETS CubicInterpolation EXPORT CubicInterpolationTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/CubicInterpolation
    )

export(EXPORT CubicInterpolationTargets
    FILE "${CMAKE_CURRENT_BINARY_DIR}/CubicInterpolationTargets.cmake")

install(EXPORT CubicInterpolationTargets
    FILE CubicInterpolationTargets.cmake
    NAMESPACE CubicInterpolation::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/CubicInterpolation
    )

include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/CubicInterpolationConfig.cmake
    INSTALL_DESTINATION lib/cmake
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
    )

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/CubicInterpolationConfig.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/CubicInterpolation
    )

install(FILES "${PROJECT_BINARY_DIR}/conanbuildinfo.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/CubicInterpolation)
