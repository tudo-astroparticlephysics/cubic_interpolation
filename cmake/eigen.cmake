include(FetchContent)

message(CHECK_START "Searching Eigen3")
find_package (Eigen3 3.3 QUIET NO_MODULE)
if(Eigen3_FOUND)
    message(STATUS "  found eigen in: ${Eigen3_DIR}")
    message(CHECK_PASS "done.")
    return()
endif()

message(STATUS "  No Eigen3 found.")
message(STATUS "  Fetching Eigen3. This may take a while.")
FetchContent_Declare(
    eigen3
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
    )

set(BUILD_TESTING OFF)
FetchContent_MakeAvailable(eigen3)

set(Eigen3_DIR ${eigen3_BINARY_DIR} CACHE PATH "Preferred Eigen3 installation prefix")

find_package (Eigen3 3.3 REQUIRED QUIET NO_MODULE)
message(STATUS "  found eigen in: ${Eigen3_DIR}")
message(CHECK_PASS "done.")
