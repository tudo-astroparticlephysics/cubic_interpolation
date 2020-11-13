include(FetchContent)
message(CHECK_START "Searching boost")

find_package(Boost 1.58.0 COMPONENTS ${BOOST_LIBS})
if(Boost_FOUND)
    message(STATUS "  found includes in: ${Boost_INCLUDE_DIRS}")
    message(STATUS "  found includes in: ${Boost_LIBRARY_DIRS}")
    message(CHECK_PASS "done.")
    return()
endif()

message(STATUS "  No boost found.")
FetchContent_Declare(
    boost
    URL https://dl.bintray.com/boostorg/release/1.74.0/source/boost_1_74_0.tar.gz
    URL_HASH SHA256=afff36d392885120bcac079148c177d1f6f7730ec3d47233aa51b0afa4db94a5
    )
FetchContent_GetProperties(boost)
if(NOT boost_POPULATED)
    message(STATUS "  Fetching Boost. This may take a while.")
    FetchContent_Populate(boost)
endif()

message(STATUS "CMAKE_INSTALL_PREFIX =${CMAKE_INSTALL_PREFIX}")
execute_process(
    COMMAND ${boost_SOURCE_DIR}/bootstrap.sh ${CMAKE_CXX_COMPILER_NAME}
    WORKING_DIRECTORY ${boost_SOURCE_DIR}
    )

set(CMAKE_CXX_COMPILER_NAME:STRING="NotFound")
if (CMAKE_CXX_COMPILER_ID MATCHES "AppleClang")
    set(CMAKE_CXX_COMPILER_NAME "clang")
    string(APPEND CMAKE_CXX_FLAGS " --gcc-toolchain")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_COMPILER_NAME "gcc")
else ()
    message( SEND_ERROR "Compiler is not supported." )
endif()


set(BUILD_FLAGS "")

foreach(lib IN LISTS BOOST_LIBS)
    list(APPEND BUILD_FLAGS "--with-${lib}")
endforeach()

execute_process(
    COMMAND ./b2 install ${BUILD_FLAGS} --prefix=${boost_BINARY_DIR}
    WORKING_DIRECTORY ${boost_SOURCE_DIR}
    )

set(BOOST_ROOT ${boost_BINARY_DIR} CACHE
    PATH "Preferred Boost installation prefix")
set(BOOST_INCLUDEDIR ${boost_BINARY_DIR}/include CACHE
    PATH "Preferred Boost include directory")
set(BOOST_LIBRARYDIR ${boost_BINARY_DIR}/lib CACHE
    PATH "Preferred Boost library directory ")
