project(${CMAKE_PROJECT_NAME})

set(sources
        main.cpp
        EdgeList.cpp
        )

set(executable "unittests")

add_executable(${executable} "")

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
        )

target_include_directories(${executable}
        PRIVATE
        ./
        PUBLIC
        "${GTEST_INCLUDE_DIRS}"
        )

add_test(NAME ${executable}
        COMMAND ${executable}
        WORKING_DIRECTORY ${binDir}
        )

