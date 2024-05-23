#include <mpi.h>
#include <unistd.h>
#include <random>

#include "Matrix.hpp"
#include "BiCGstab.hpp"

using std::vector;

template<typename T>
void gen_matrix(Matrix<T>& A){
    std::default_random_engine eng;
    std::uniform_real_distribution<double> distr(-10, 10);
    int n = A.size().first, m = A.size().second;
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < m; ++j){
            A(i, j) = distr(eng);
        }
    }
    for(int i = 0; i < n; ++i){
        A(i, i) = fabs(A(i, i));
        for(int j = 0; j < m; ++j){
            if(i != j){
                A(i, i) += sqrt(fabs(A(i, j)));
            }
        }
    }
}


int main(int argc, char**argv){
    int mpi_thread_prov;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &mpi_thread_prov);
    int proc_id, n_proc;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);

    int N = atoi(argv[1]);
    int k = N * n_proc;
    Matrix<double> A(proc_id, n_proc, MPI_COMM_WORLD), A_chunk(proc_id, n_proc, MPI_COMM_WORLD);
    std::vector<double> b, b_chunk;
    if(proc_id == 0){
        int n = N * n_proc;
        A.resize(n, k);
        gen_matrix(A);
        b.resize(k);
        for(int i = 0; i < k; ++i){
            b[i] = sin(i);
        }
        if(n < 10) {
            print_matrix(A);
            std::cout << "-----------------\n";
        }
    }
    A_chunk.resize(N, k);
    b_chunk.resize(N);
    MPI_Scatter(A.data(), N * k, MPI_DOUBLE, A_chunk.data(), N * k, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(b.data(), N, MPI_DOUBLE, b_chunk.data(), N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //cblas_dgemv(CblasRowMajor, CblasNoTrans, N, N, 1, &(A_chunk(0, N)), k, b_chunk.data(), 1, 0, x_chunk.data(), 1);

    vector<double> x = BiCGstab<Matrix<double>>(A_chunk, b_chunk, 1e-8);

    MPI_Barrier(MPI_COMM_WORLD);
    vector<double> r = b_chunk - A_chunk.matvec(x);
    std::cout << "r = " << norm(r, proc_id, n_proc, MPI_COMM_WORLD) << std::endl;
    MPI_Finalize();
}