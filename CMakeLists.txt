cmake_minimum_required(VERSION 3.20)
project(FDM)


if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "release")
endif()

set(cmake_cxx_flags_release "-Ofast -flto -march=native")

find_package(inmost REQUIRED)

link_directories(${INMOST_LIBRARY_DIRS})
include_directories(${INMOST_INCLUDE_DIRS})
add_definitions(${INMOST_DEFINITIONS})


add_executable(fdm main.cpp
        functions.cpp)
target_link_libraries(fdm ${INMOST_LIBRARIES})