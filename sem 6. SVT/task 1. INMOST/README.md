# Compilation
## Requriements:
1. INMOST package

## Build
```console
mkdir build && cd build
cmake .. && make
```
## Usage
```console
./svt ../data/A_g2_1.mtx ../data/rhs_g2_1.mtx 1e-7
```
3rd argument is a drop tolerance for Solver