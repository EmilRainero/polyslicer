cmake_minimum_required(VERSION 3.10)
project(polyslicer)

set(CMAKE_CXX_STANDARD 11)

set( sources
        main.cpp
        )

set( headers
        Vector3.h
        LineSegment.h
        Triangle.h
        TriangleMesh.h
        Timer.h
        )

add_executable(polyslicer ${sources})