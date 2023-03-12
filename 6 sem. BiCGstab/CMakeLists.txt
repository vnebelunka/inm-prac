cmake_minimum_required(VERSION 3.22)
project(BiCGstab)

set(CMAKE_CXX_STANDARD 20)
set(GCC_COVERAGE_COMPILE_FLAGS "-fsanitize=address")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS} -lcblas")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

find_package(BLAS REQUIRED)


add_executable(BiCGstab main.cpp Matrix.cpp)
target_include_directories(BiCGstab PRIVATE include)