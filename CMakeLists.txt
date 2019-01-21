cmake_minimum_required(VERSION 3.10)
project(polyslicer)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")

set(CMAKE_CXX_STANDARD 11)

set(sources
        main.cpp
        PolygonSlicer.cpp
        Image.cpp
        EdgeTable.cpp
        Point.cpp
        Polygon.cpp
        Node.cpp
        LodePNGImage.cpp
        LodePNG.cpp
        )

set( headers
        Vector3.h
        LineSegment.h
        Triangle.h
        TriangleMesh.h
        Timer.h
        )

FIND_PACKAGE (Threads)

#add_executable(polyslicer ${sources})

set( executable "polyslicer" )
add_executable( ${executable} "")

target_link_libraries( ${executable}  ${CMAKE_THREAD_LIBS_INIT} )

set_target_properties(${outputLib} PROPERTIES
        CXX_STANDARD 11
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS ON
        )

target_sources( ${executable}
        PRIVATE
        ${sources}
        ${headers}
        )

#include(GoogleTest)
#
#add_subdirectory("unittests")

cmake_minimum_required(VERSION 2.6)

# Locate GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_subdirectory("unittests")

# Link runTests with what we want to test and the GTest and pthread library
#add_executable(executeTests sqrt_test.cpp)
#target_link_libraries(executeTests ${GTEST_LIBRARIES} pthread)
