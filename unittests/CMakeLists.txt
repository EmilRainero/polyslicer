project(${CMAKE_PROJECT_NAME})

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")

set(sources
        main.cpp
        EdgeListTest.cpp

        ../Image.cpp
        ../EdgeTable.cpp
        )

set(executable "unittests")

add_executable(${executable} ${sources})

set_target_properties(${outputLib} PROPERTIES
        CXX_STANDARD 11
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS ON
        )

target_sources(${executable}
        PRIVATE
        ${sources}
        )

target_link_libraries(${executable}
        PUBLIC
        "${GTEST_LIBRARIES}"
        pthread
        )

target_include_directories(${executable}
        PUBLIC
        ../
        "${GTEST_INCLUDE_DIRS}"
        )

