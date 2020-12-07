set(CUBIC_INTERPOLATION_SRC
    Axis.cxx
    BicubicSplines.cxx
    CubicSplines.cxx
    InterpolantBuilder.cxx
    Interpolant.cxx
    )

target_sources(CubicInterpolation PRIVATE ${CUBIC_INTERPOLATION_SRC})

target_include_directories(CubicInterpolation PUBLIC
    ${Boost_INCLUDE_DIRS}
    )

target_link_libraries(CubicInterpolation
    PUBLIC
        Eigen3::Eigen
    PRIVATE
        Boost::serialization
        Boost::filesystem
    )
