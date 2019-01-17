cmake_minimum_required(VERSION 3.10)
project(polyslicer)

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")

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
