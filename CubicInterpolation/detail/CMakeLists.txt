set(CUBIC_INTERPOLATION_SRC
    BicubicSplines.cxx
    CubicSplines.cxx
    InterpolantBuilder.cxx
    Interpolant.cxx
    )

set(CUBIC_INTERPOLATION_HEADERS
    Interpolant.hpp
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

install(FILES ${CUBIC_INTERPOLATION_HEADERS}
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/CubicInterpolation/detail/")
