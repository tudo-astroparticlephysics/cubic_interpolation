configure_file(version.h.in version.h @ONLY)

set (CubicInterpolation_header
    ${CMAKE_CURRENT_BINARY_DIR}/version.h
    CubicInterpolation/Axis.h
    CubicInterpolation/BicubicSplines.h
    CubicInterpolation/CMakeLists.txt
    CubicInterpolation/CubicSplines.h
    CubicInterpolation/FindParameter.hpp
    CubicInterpolation/Interpolant.h
    CubicInterpolation/Interpolant.hpp
    CubicInterpolation/InterpolantBuilder.h
    )

set_target_properties(CubicInterpolation PROPERTIES
    PUBLIC_HEADER "${CubicInterpolation_header}")
