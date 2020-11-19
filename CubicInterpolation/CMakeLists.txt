add_library(CubicInterpolation INTERFACE)

target_include_directories(CubicInterpolation INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    ${Boost_INCLUDE_DIRS}
    )

target_link_libraries (
    CubicInterpolation
    INTERFACE
    Eigen3::Eigen
    Boost::serialization
    Boost::filesystem
    )


install(TARGETS CubicInterpolation
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
