#include "Matrix.hpp"
#include <iostream>
#include <cblas.h>
#include <unistd.h>
#include <mpi.h>

double dot( std::vector<double> &a_chunk, std::vector<double> & b_chunk, int proc_id, int n_proc){
    double dot_chunk = cblas_ddot(a_chunk.size(), a_chunk.data(), 1, b_chunk.data(), 1);
    double dot;
    MPI_Allreduce(&dot_chunk, &dot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return dot;
}

double norm(std::vector<double> &a_chunk, int proc_id, int n_proc){
    return sqrt(dot(a_chunk, a_chunk, proc_id, n_proc));
}


std::vector<double> matvec(Matrix<double> &A_chunk, const std::vector<double> &x_chunk, int proc_id, int n_proc){
    int k = A_chunk.size().second;
    int n = A_chunk.size().first;
    std::vector<double> x(n * n_proc), y_chunk(n);
    MPI_Allgather(x_chunk.data(), n, MPI_DOUBLE, x.data(), n, MPI_DOUBLE, MPI_COMM_WORLD);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, n, k, 1, A_chunk.data(), k, x.data(), 1, 0, y_chunk.data(), 1);
    return y_chunk;
}


void print_matrix(Matrix<double> const& A){
    size_t n = A.size().first;
    size_t m = A.size().second;
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < m; ++j){
            std::cout << A(i, j) << ' ';
        }
        std::cout << '\n';
    }
}

void print_matrix(Matrix<double> const& A, int proc_id){
    sleep(proc_id);
    size_t n = A.size().first;
    size_t m = A.size().second;
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < m; ++j){
            std::cout << A(i, j) << ' ';
        }
        std::cout << '\n';
    }
}
