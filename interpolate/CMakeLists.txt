add_library(INTERPOLATE INTERFACE)

target_include_directories(INTERPOLATE INTERFACE
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_sources(INTERPOLATE INTERFACE
    BicubicSplines.h
    Utility.h
    Axis.h
    InterpolantBuilder.h
    )

find_package(Eigen3)
find_package(Boost QUIET COMPONENTS serialization)

target_link_libraries (
  INTERPOLATE
  INTERFACE
  Eigen3::Eigen
  Boost::serialization
  )

install(TARGETS INTERPOLATE
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
