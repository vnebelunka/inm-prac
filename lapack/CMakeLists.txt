project(Cucumber_blas)

set(GCC_COVERAGE_COMPILE_FLAGS "-flto -Ofast -march=native -lcblas")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

find_package(LAPACK REQUIRED)


add_executable(lapack main.cpp Matrix.cpp Matrix.hpp)