cmake_minimum_required(VERSION 3.22)
project(BiCGstab)

set(CMAKE_CXX_STANDARD 20)


set(CMAKE_CXX_FLAGS_DEBUG_INIT "-Wall -fsanitize=address -Wextra -lcblas")
set(CMAKE_CXX_FLAGS_RELEASE_INIT "-Ofast -march=native -flto -lcblas")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lcblas")

find_package(BLAS REQUIRED)
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

add_executable(BiCGstab main.cpp Matrix.cpp)
target_link_libraries(BiCGstab OpenMP::OpenMP_CXX MPI::MPI_CXX)
target_include_directories(BiCGstab PUBLIC include)

add_executable(BiCGstab_mpi mpi.cpp Matrix.cpp)
target_include_directories(BiCGstab_mpi PUBLIC include)
target_link_libraries(BiCGstab_mpi MPI::MPI_CXX)