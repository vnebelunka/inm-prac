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


/*std::vector<double> matvec(Matrix<double> &A_chunk, const std::vector<double> &x_chunk, int proc_id, int n_proc){
    int k = A_chunk.size().second;
    int n = A_chunk.size().first;
    std::vector<double> x(n * n_proc), y_chunk(n);
    MPI_Allgather(x_chunk.data(), n, MPI_DOUBLE, x.data(), n, MPI_DOUBLE, MPI_COMM_WORLD);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, n, k, 1, A_chunk.data(), k, x.data(), 1, 0, y_chunk.data(), 1);
    return y_chunk;
}*/

std::vector<double> matvec(Matrix<double> &A_chunk, const std::vector<double> &x_chunk, int proc_id, int n_proc){
    int k = A_chunk.size().second;
    int n = A_chunk.size().first;
    std::vector<double> y_chunk(n), matvec_chunk(x_chunk);
    for(int shift = 1; shift < n_proc; ++shift){
        int send_proc_id = (proc_id + shift) % n_proc, recv_proc_id = (proc_id - shift + n_proc) % n_proc;
        std::vector<double> next_chunk(n);
        MPI_Request send, recv;
        MPI_Isend(x_chunk.data(), n, MPI_DOUBLE, send_proc_id, shift, MPI_COMM_WORLD, &send);
        MPI_Irecv(next_chunk.data(), n, MPI_DOUBLE, recv_proc_id, shift, MPI_COMM_WORLD, &recv);

        int block_start_col = (proc_id + shift - 1) % n_proc;
        cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1, &A_chunk(0, block_start_col * n), k, matvec_chunk.data(), 1, 1, y_chunk.data(), 1);

        MPI_Status status;
        MPI_Wait(&send, &status);
        if(!status.MPI_ERROR){
            throw std::runtime_error("send failed");
        }
        MPI_Wait(&recv, &status);
        if(!status.MPI_ERROR){
            throw std::runtime_error("recv failed");
        }
        matvec_chunk = next_chunk;
    }
    int block_start_col = (proc_id - 1 + n_proc) % n_proc;
    cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1, &A_chunk(0, block_start_col * n), k, matvec_chunk.data(), 1, 1, y_chunk.data(), 1);
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
