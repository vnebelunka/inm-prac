
#include "Matrix.hpp"
#include <iostream>

void print_matrix(Matrix<double> const& A){
    size_t n = A.size().first;
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            printf("%4lf ", A(i, j));
        }
        putchar('\n');
    }
    putchar('\n');
}