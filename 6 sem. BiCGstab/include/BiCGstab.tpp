#pragma once

#include <vector>
#include "Matrix.hpp"
#include "BiCGstab.hpp"

using std::vector;

template<typename matrix>
vector<double> BiCGstab(const matrix &A_chunk, vector<double>& b_chunk, double rel_tol){
    int n = b_chunk.size();
    int n_proc = A_chunk.get_n_proc(), proc_id = A_chunk.get_proc_id();
    MPI_Comm comm = A_chunk.get_comm();
    rel_tol *= norm(b_chunk, proc_id, n_proc, comm);
    vector<double> x(n), r_k = b_chunk - A_chunk.matvec(x), rconj0 = r_k, p = r_k;
    int iter = 0;
    for(;iter < 1e4; ++iter){
        vector<double> Ap_k = A_chunk.matvec(p);
        double alpha_k = dot(r_k, rconj0, proc_id, n_proc, comm) / dot(Ap_k, rconj0, proc_id, n_proc, comm);
        vector<double> s_k = r_k - alpha_k * Ap_k;
        if(norm(s_k, proc_id, n_proc, comm) < rel_tol){
            x = x + alpha_k * p;
            break;
        }
        vector<double> As_k = A_chunk.matvec(s_k);
        double w_k = dot(As_k, s_k, proc_id, n_proc, comm) / dot(As_k, As_k, proc_id, n_proc, comm);
        x = x + alpha_k * p + w_k * s_k;
        double r_krconj_0 = dot(r_k, rconj0, proc_id, n_proc, comm);
        // orthogonality check
        if(iter % 10 != 1) {
            r_k = s_k - w_k * As_k;
        } else {
            r_k = b_chunk - A_chunk.matvec(x);
            std::cout << "r norm = " << norm(r_k, proc_id, n_proc, comm) << " on iteration "  << iter << std::endl;
        }
        if(norm(r_k, proc_id, n_proc, comm) < rel_tol){
            break;
        }
        double b_k = dot(r_k, rconj0, proc_id, n_proc, comm) / r_krconj_0;
        p = r_k + b_k * p - w_k * Ap_k;
        // restart
        if(fabs(dot(r_k, rconj0, proc_id, n_proc, comm)) < 1e-8){
            rconj0 = r_k;
            p = r_k;
        }
    }
    std::cout << "iters: " << iter << std::endl;
    return x;
}