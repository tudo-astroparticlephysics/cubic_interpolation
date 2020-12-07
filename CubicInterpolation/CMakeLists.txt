add_subdirectory(details)

set(CUBIC_INTERPOLATION_HEADERS
    Axis.h
    BicubicSplines.h
    CubicSplines.h
    Interpolant.h
    InterpolantBuilder.h
    )

target_include_directories(CubicInterpolation
    PRIVATE
        ${Boost_INCLUDE_DIRS}
    )

install(FILES ${CUBIC_INTERPOLATION_HEADERS}
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/CubicInterpolation")
