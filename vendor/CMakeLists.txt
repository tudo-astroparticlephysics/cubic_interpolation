message(STATUS "Finding Eigen3")
find_package(Eigen3 3.3 NO_MODULE)

if (Eigen3_FOUND)
    message(STATUS "found in: ${EIGEN3_ROOT_DIR}")
else()
    message(STATUS "No Eigen3 installation found.")
    include(FetchContent)
    FetchContent_Declare(
      Eigen3
      GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
      GIT_SHALLOW TRUE
      GIT_PROGRESS TRUE
    )

    FetchContent_GetProperties(Eigen3)
    if(NOT eigen3_POPULATED)
      message(STATUS "Fetching Eigen3 master.")
      FetchContent_Populate(eigen3)

      add_subdirectory(${eigen3_SOURCE_DIR} ${eigen3_BINARY_DIR})
    endif()
endif ()
