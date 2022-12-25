cmake_minimum_required(VERSION 3.22)
project(Cucumber)

set(CMAKE_CXX_STANDARD 20)

set(GCC_COVERAGE_COMPILE_FLAGS "-O2 -fsanitize=address")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

add_subdirectory(blas)
add_subdirectory(fortran)
add_subdirectory(lapack)

add_executable(Cucumber main.cpp progressbar.h Matrix.cpp Matrix.hpp blas/main.cpp)
