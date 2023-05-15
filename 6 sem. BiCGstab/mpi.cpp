#include "Matrix.hpp"
#include <mpi.h>
#include <unistd.h>
#include <random>
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



vector<double> BiCGstab(Matrix<double> &A_chunk, vector<double>& b_chunk, double eps, int proc_id, int n_proc){
    int n = b_chunk.size();
    vector<double> x(n), r_k = b_chunk - matvec(A_chunk, x, proc_id, n_proc), rconj0 = r_k, p = r_k;
    int iter = 0;
    for(;iter < 1e4; ++iter){
        vector<double> Ap_k = matvec(A_chunk, p, proc_id, n_proc);
        double alpha_k = dot(r_k, rconj0, proc_id, n_proc) / dot(Ap_k, rconj0, proc_id, n_proc);
        vector<double> s_k = r_k - alpha_k * Ap_k;
        if(norm(s_k, proc_id, n_proc) < eps){
            x = x + alpha_k * p;
            break;
        }
        vector<double> As_k = matvec(A_chunk, s_k, proc_id, n_proc);
        double w_k = dot(As_k, s_k, proc_id, n_proc) / dot(As_k, As_k, proc_id, n_proc);
        x = x + alpha_k * p + w_k * s_k;
        double r_krconj_0 = dot(r_k, rconj0, proc_id, n_proc);
        // orthogonality check
        if(iter % 1000 != 1) {
            r_k = s_k - w_k * As_k;
        } else {
            r_k = b_chunk - matvec(A_chunk, x, proc_id, n_proc);
            std::cout << norm(r_k, proc_id, n_proc) << std::endl;
        }
        if(norm(r_k, proc_id, n_proc) < eps){
            break;
        }
        double b_k = dot(r_k, rconj0, proc_id, n_proc) / r_krconj_0;
        p = r_k + b_k * p - w_k * Ap_k;
        // restart
        if(fabs(dot(r_k, rconj0, proc_id, n_proc)) < 1e-8){
            rconj0 = r_k;
            p = r_k;
        }
    }
    std::cout << "iters: " << iter << std::endl;
    return x;
}


int main(int argc, char**argv){
    int mpi_thread_prov;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &mpi_thread_prov);
    int proc_id, n_proc;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);

    int N = atoi(argv[1]);
    int k = N * n_proc;
    Matrix<double> A, A_chunk;
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

    vector<double> x = BiCGstab(A_chunk, b_chunk, 1e-5, proc_id, n_proc);
    MPI_Barrier(MPI_COMM_WORLD);

    vector<double> r = b_chunk - matvec(A_chunk, x, proc_id, n_proc);
    std::cout << "r = " << norm(r, proc_id, n_proc) << std::endl;
    MPI_Finalize();
}