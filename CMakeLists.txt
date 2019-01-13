cmake_minimum_required(VERSION 3.10)
project(polyslicer)

set(CMAKE_CXX_STANDARD 11)

set( sources
        main.cpp
        PolygonSlicer.cpp
        )

set( headers
        Vector3.h
        LineSegment.h
        Triangle.h
        TriangleMesh.h
        Timer.h
        )


#add_executable(polyslicer ${sources})

set( executable "polyslicer" )
add_executable( ${executable} "")

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
