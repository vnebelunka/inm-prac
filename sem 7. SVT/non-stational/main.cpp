//
// Created by nebil on 22.10.23.
//
#include <cstdio>
#include <vector>
#include <stdexcept>
#include <iostream>

#include "inmost.h"
#include "functions.hpp"

using namespace INMOST;


void save_solution(size_t n, const Sparse::Vector &sol, std::string vtk_fname){
    size_t N = n * n;
    double h = 1. / n;
    FILE* fd = fopen(vtk_fname.c_str(), "w");
    if(!fd){
        throw std::runtime_error("can't open file " + vtk_fname);
    }
    fprintf(fd, "# vtk DataFile Version 2.0\n");
    fprintf(fd, "solution by FDM\n");
    fprintf(fd, "ASCII\nDATASET UNSTRUCTURED_GRID\n");
    fprintf(fd, "POINTS %lu double\n", N);
    for(int i = 0; i * i < N; ++i){
        for(int j = 0; j * j < N; ++j) {
            fprintf(fd, "%lf %lf %lf\n", i * h, j * h, 0.0);
        }
    }
    fprintf(fd, "POINT_DATA %lu\n",  N);
    fprintf(fd, "SCALARS U double 1\n");
    fprintf(fd, "LOOKUP_TABLE default\n");
    for(int i = 0; i < N; ++i){
        fprintf(fd, "%lf\n", sol[i]);
    }
}

void calc_matrix(Sparse::Matrix& A, const size_t n, const double h){
    size_t N = n * n;
    A.SetInterval(0, N);
    const double h2 = h * h;
#pragma omp parallel for
    for(int k = 0; k < N; ++k){
        int i = k / n, j = k % n;
        double x = h * i, y = h * j;
        if(i == 0 || j == 0 || i == n - 1 || j == n - 1){ // boundary points
            A[k][k] = 1;
        } else {
            A[k][k] = 2 * (dx + dy) / h2;
            if(Problem_type == IMPLICIT){
                A[k][k] += 1. / dt;
            }
            A[k][k + n] = -dx / h2;
            A[k][k - n] = -dx / h2;
            A[k][k - 1] = -dy / h2;
            A[k][k + 1] = -dy / h2;
            A[k][k + n + 1] = -dxy / (2 * h2);
            A[k][k - n - 1] = -dxy / (2 * h2);
            A[k][k - n + 1] = dxy / (2 * h2);
            A[k][k + n - 1] = dxy / (2 * h2);
        }
    }
}


void apply_iteration(Solver& S, const size_t n, double h, Sparse::Vector &u, int t_step, Sparse::Matrix& A){
    size_t N = n * n;
    if(Problem_type == IMPLICIT){
        double t = t_step * dt;
        for(int k = 0; k < N; ++k){
            int i = k / n, j = k % n;
            double x = h * i, y = h * j;
            if(i == 0 || j == 0 || i == n - 1 || j == n - 1){ // boundary points
                u[k] = phi(x, y, t);
            } else {
                u[k] /= dt;
                u[k] += f(x, y, t);
            }
        }
        bool solved = S.Solve(u, u);
        if(!solved){
            std::cout << "Linear solver failure!" << std::endl;
            std::cout << "Reason: " << S.ReturnReason() << std::endl;
        }
    } else {
        double t = (t_step - 1) * dt;
        Sparse::Vector temp;
        temp.SetInterval(0, N);
        A.MatVec(1.0, u, 0.0, temp);
        for(int k = 0; k < N; ++k){
            int i = k / n, j = k % n;
            double x = h * i, y = h * j;
            u[k] = u[k] + (f(x, y, t) - temp[k]) * dt;
        }
    }
    save_solution(n, u, std::string("res") + std::to_string(t_step) + std::string(".vtk"));
}

int main(int argc, char**argv){
    Solver::Initialize(&argc, &argv);
    int n = strtol(argv[1], nullptr, 10), N = n * n;
    double h = 1. / n;
    std::cout << "h: " << h << " n: " << n << std::endl;

    std::cout << "dx: " << dx << " dy: " << dy << std::endl;

    Sparse::Matrix A;
    calc_matrix(A, n, h);
    Sparse::Vector u;
   
    u.SetInterval(0, N);
    for(int k = 0; k < N; ++k){
        int i = k / n, j = k % n;
        double x = h * i, y = h * j;
        u[k] = u0(x, y);
    }
    save_solution(n, u, "res0.vtk");

    Solver S(Solver::INNER_ILU2);
    S.SetParameter("absolute_tolerance", "1e-12");
    S.SetParameter("verbosity", "2");
    S.SetParameter("relative_tolerance", "1e-7");
    S.SetMatrix(A);
    
    for(int t_step = 1; t_step <= T; ++t_step){
        apply_iteration(S, n, h, u, t_step, A);
    }

    Solver::Finalize();
}