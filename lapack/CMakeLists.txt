cmake_minimum_required(VERSION 3.22)
project(Cucumber)

enable_language(Fortran)

set(GCC_COVERAGE_COMPILE_FLAGS "-flto -Ofast -march=native -llapack")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
set(CMAKE_Fortran_FLAGS "${GCC_COVERAGE_COMPILE_FLAGS}")

find_package(LAPACK REQUIRED)


add_executable(lapack main.cpp Matrix.cpp Matrix.hpp)
add_executable(lapackf lapack.f90)