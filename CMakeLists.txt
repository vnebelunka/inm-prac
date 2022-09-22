cmake_minimum_required(VERSION 3.22)
project(Cucumber)

set(CMAKE_CXX_STANDARD 20)

set(GCC_COVERAGE_COMPILE_FLAGS "-O2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")


add_executable(Cucumber main.cpp progressbar.h Matrix.cpp Matrix.hpp)
