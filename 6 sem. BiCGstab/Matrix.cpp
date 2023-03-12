#include "Matrix.hpp"
#include <iostream>
#include <cblas.h>

template<>
double dot(const std::vector<double>& a, const std::vector<double>& b){
    return cblas_ddot(a.size(), a.data(), 1, b.data(), 1);
}

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