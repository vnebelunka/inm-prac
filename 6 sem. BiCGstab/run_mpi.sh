#!/bin/bash

mkdir -p build
cd build || { echo "Can not create build directory"; exit 1; }


export OMP_NUM_THREADS=4
cmake ..
make -j BiCGstab_mpi && mpirun -np 4 ./BiCGstab_mpi 3000
