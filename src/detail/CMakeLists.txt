target_sources(CubicInterpolation PRIVATE
    Axis.cxx
    BicubicSplines.cxx
    BicubicSplines_imp.cxx
    CubicSplines.cxx
    CubicSplines_imp.cxx
    InterpolantBuilder.cxx
    Interpolant.cxx
    )

find_package(Boost REQUIRED serialization filesystem)
find_package(Eigen3 REQUIRED NO_MODULE)

target_link_libraries(CubicInterpolation
    PRIVATE
    Eigen3::Eigen
    Boost::serialization
    Boost::filesystem
    )
