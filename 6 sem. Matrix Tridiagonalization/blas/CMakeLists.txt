project(Cucumber)

enable_language(Fortran)

set(GCC_COVERAGE_COMPILE_FLAGS "-flto -Ofast -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS} -lcblas")
set(CMAKE_Fortran_FLAGS "${GCC_COVERAGE_COMPILE_FLAGS} -lblas")

find_package(BLAS REQUIRED)

add_executable(blas main.cpp Matrix.cpp Matrix.hpp)
add_executable(blasf blas.f90)