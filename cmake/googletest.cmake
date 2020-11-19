include(FetchContent)

message(CHECK_START "Searching googletest")
find_package(GTest COMPONENTS GTest::GTest)
if(GTest_FOUND)
    message(STATUS "  found includes in: ${GTEST_INCLUDE_DIRS}")
    message(STATUS "  found lib in: ${GTEST_LIBRARIES}")
    message(CHECK_PASS "done.")
    return()
endif()

message(STATUS "  No GTest found.")
message(STATUS "  No googletest found.")
message(STATUS "  Fetching googletest. This may take a while.")
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
    )
FetchContent_MakeAvailable(googletest)
message(CHECK_PASS "done.")
