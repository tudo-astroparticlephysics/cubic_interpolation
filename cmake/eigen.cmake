include(FetchContent)

message(CHECK_START "Searching Eigen3")
find_package (Eigen3 3.3 QUIET NO_MODULE)
if (NOT Eigen3_FOUND)
    message(STATUS "  No Eigen3 found.")
    message(STATUS "  Fetching Eigen3. This may take a while.")
    FetchContent_Declare(
        Eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_SHALLOW TRUE
        GIT_PROGRESS TRUE
        )
    FetchContent_MakeAvailable(Eigen3)
endif()
message(STATUS "  found includes in: ${EIGEN3_INCLUDE_DIRS}")
message(CHECK_PASS "done.")
