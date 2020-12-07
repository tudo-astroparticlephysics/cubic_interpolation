include(FetchContent)

message(CHECK_START "Searching boost")
find_package(Boost 1.58.0 QUIET COMPONENTS ${BOOST_LIBS})
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

execute_process(
    COMMAND ${boost_SOURCE_DIR}/bootstrap.sh ${CMAKE_CXX_COMPILER_NAME}
    WORKING_DIRECTORY ${boost_SOURCE_DIR}
    )

set(CMAKE_CXX_COMPILER_NAME:STRING="NotFound")
if (CMAKE_CXX_COMPILER_ID MATCHES "AppleClang")
    set(CMAKE_CXX_COMPILER_NAME "clang")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_COMPILER_NAME "gcc")
else ()
    message( SEND_ERROR "Compiler is not supported." )
endif()


set(BUILD_FLAGS "")

foreach(lib IN LISTS BOOST_LIBS)
    list(APPEND BUILD_FLAGS "--with-${lib}")
endforeach()

if(Boost_USE_STATIC_RUNTIME)
    list(APPEND BUILD_FLAGS "runtime-link=static")
endif()

if(Boost_USE_DEBUG_RUNTIME)
    list(APPEND BUILD_FLAGS "runtime-debugging=on")
endif()

if(Boost_USE_DEBUG_PYTHON)
    list(APPEND BUILD_FLAGS "python-debugging=on")
endif()

if(Boost_USE_DEBUG_LIBS)
    list(APPEND BUILD_FLAGS "variant=debug")
endif()

if(Boost_USE_STLPORT)
    list(APPEND BUILD_FLAGS "stdlib=stlport")
else()
    list(APPEND BUILD_FLAGS "stdlib=native")
endif()

foreach(flag IN LISTS BOOST_FLAGS)
    message(STATUS "BOOST_FLAG: ${flag}")
endforeach()

execute_process(
    COMMAND ./b2 install ${BUILD_FLAGS} --prefix=${CMAKE_INSTALL_PREFIX}
    WORKING_DIRECTORY ${boost_SOURCE_DIR}
    )

set(BOOST_ROOT ${CMAKE_INSTALL_PREFIX} CACHE PATH "Preferred Boost installation prefix")
set(BOOST_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR} CACHE PATH "Preferred Boost include directory")
set(BOOST_LIBRARYDIR ${CMAKE_INSTALL_LIBDIR} PATH "Preferred Boost library directory ")

find_package(Boost 1.58.0 REQUIRED QUIET COMPONENTS ${BOOST_LIBS})

message(CHECK_PASS "done.")
