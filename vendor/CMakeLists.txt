# SET(FETCHCONTENT_QUIET FALSE)
include(FetchContent)

message(CHECK_START "Searching Eigen3")
find_package (Eigen3 3.3 NO_MODULE)

if (Eigen3_FOUND)
    message(STATUS "  found includes in: ${EIGEN3_INCLUDE_DIRS}")
    message(CHECK_PASS "done.")
else()
    message(CHECK_FAIL "No Eigen3 found.")
    message(STATUS "  Fetching Eigen3. This may take a while.")
    FetchContent_Declare(
        Eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_SHALLOW TRUE
        GIT_PROGRESS TRUE
        )
    FetchContent_MakeAvailable(Eigen3)
endif()

set(BOOST_LIBS serialization)

message(CHECK_START "Searching boost")
find_package(Boost 1.58.0 QUIET COMPONENTS ${BOOST_LIBS})
if(Boost_FOUND)
    message(STATUS "  found includes in: ${boost_INCLUDE_DIRS}")
    message(STATUS "  found includes in: ${boost_LIBRARY_DIRS}")
    message(CHECK_PASS "done.")
else()
    message(CHECK_FAIL "No boost found.")
    FetchContent_Declare(
        boost
        URL https://dl.bintray.com/boostorg/release/1.74.0/source/boost_1_74_0.tar.gz
        URL_HASH SHA256=afff36d392885120bcac079148c177d1f6f7730ec3d47233aa51b0afa4db94a5
        )
    FetchContent_GetProperties(boost)
    if(NOT boost_POPULATED)
        message(STATUS "  Fetching Boost. This may take a while.")
        FetchContent_Populate(boost)

        option(ENABLE_BOOST_BUILD "Enable reconfiguring and rerunning the boost build" ON)
        if(ENABLE_BOOST_BUILD)

            # This file comes from the following location:
            #   https://github.com/pfultz2/cget/blob/master/cget/cmake/boost.cmake
            configure_file(boost.cmake
                ${boost_SOURCE_DIR}/CMakeLists.txt
                COPYONLY
                )

            unset(generatorArgs)
            set(cacheArgs
                "-DCMAKE_INSTALL_PREFIX:STRING=${boost_BINARY_DIR}/install"
                "-DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=YES"
                "-DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}"
                )
            if(CMAKE_TOOLCHAIN_FILE)
                list(APPEND cacheArgs "-DCMAKE_TOOLCHAIN_FILE:FILEPATH=${CMAKE_TOOLCHAIN_FILE}")
            else()
                list(APPEND cacheArgs "-DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}"
                    "-DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}"
                    )
            endif()

            get_property(isMulti GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
            if(NOT isMulti)
                list(APPEND cacheArgs "-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}")
            endif()

            if(CMAKE_GENERATOR_PLATFORM)
                list(APPEND generatorArgs
                    --build-generator-platform "${CMAKE_GENERATOR_PLATFORM}"
                    )
            endif()
            if(CMAKE_GENERATOR_TOOLSET)
                list(APPEND generatorArgs
                    --build-generator-toolset  "${CMAKE_GENERATOR_TOOLSET}"
                    )
            endif()

            # Assume parent dir has set BOOST_LIBS to a list of the boost modules
            # we want build and made available
            foreach(lib IN LISTS BOOST_LIBS)
                list(APPEND cacheArgs "-DBOOST_WITH_${lib}:STRING=--with-${lib}")
            endforeach()

            message(STATUS "Configuring and building boost immediately")
            execute_process(
                COMMAND ${CMAKE_CTEST_COMMAND}
                --build-and-test  ${boost_SOURCE_DIR} ${boost_BINARY_DIR}
                --build-generator ${CMAKE_GENERATOR} ${generatorArgs}
                --build-target    install
                --build-noclean
                --build-options   ${cacheArgs}
                WORKING_DIRECTORY ${boost_SOURCE_DIR}
                OUTPUT_FILE       ${boost_BINARY_DIR}/build_output.log
                ERROR_FILE        ${boost_BINARY_DIR}/build_output.log
                RESULT_VARIABLE   result
                )
            message(STATUS "boost build complete")
            if(result)
                message(FATAL_ERROR "Failed boost build, see build log at:\n"
                    "    ${boost_BINARY_DIR}/build_output.log")
            endif()
        endif()
        set(CMAKE_PREFIX_PATH ${boost_BINARY_DIR}/install)
    endif()
endif()
find_package(Boost 1.74.0 COMPONENTS ${BOOST_LIBS})
message(STATUS "done.")
