message(CHECK_START "Searching Eigen3")
find_package(Eigen3 3.3 NO_MODULE)

if (Eigen3_FOUND)
    message(STATUS "  found includes in: ${EIGEN3_INCLUDE_DIRS}")
    message(CHECK_PASS "done.")
else()
    message(CHECK_FAIL "No Eigen3 installation found.")
    include(FetchContent)
    FetchContent_Declare(
      Eigen3
      GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
      GIT_SHALLOW TRUE
      GIT_PROGRESS TRUE
    )

    FetchContent_MakeAvailable(Eigen3)
endif()


message(CHECK_START "Searching Boost")
find_package(Boost 1.56 QUIET COMPONENTS
    serialization)
if(Boost_FOUND)
    message(STATUS "  found includes in: ${Boost_INCLUDE_DIRS}")
    message(STATUS "  found includes in: ${Boost_LIBRARY_DIRS}")
    message(CHECK_PASS "done.")
else()
    message(STATUS "No BOOST installation found.")
endif()

message(STATUS "done.")
