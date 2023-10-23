//
// Created by nebil on 22.10.23.
//
#include <cstdio>
#include <vector>
#include <stdexcept>
#include <iostream>

#include "inmost.h"

using namespace INMOST;

extern double f(double x, double y);
extern double g_bond(double x, double y);

extern double u(double x, double y, double dx = 1., double dy = 1.);


void save_solution(size_t n, const std::vector<double>& sol, const std::vector<double>& exact_sol, std::string vtk_fname){
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
    fputc('\n', fd);
    if(!exact_sol.empty()) {
        fprintf(fd, "SCALARS U_true double 1\n");
        fprintf(fd, "LOOKUP_TABLE default\n");
        for (int i = 0; i < N; ++i) {
            fprintf(fd, "%lf\n", exact_sol[i]);
        }
        fputc('\n', fd);
    }
}


int main(int argc, char**argv){
    Solver::Initialize(&argc, &argv);
    int n = strtol(argv[1], nullptr, 10), N = n * n;
    double h = 1. / n;
    std::cout << "h: " << h << " n: " << n << std::endl;
    double dx = 1.0, dy = 1.0;
    if(argc > 3){
        dx = strtod(argv[3], nullptr);
    }
    if(argc > 4){
        dy = strtod(argv[4], nullptr);
    }
    std::cout << "dx: " << dx << " dy: " << dy << std::endl;

    Sparse::Matrix A;
    Sparse::Vector b;
    Sparse::Vector sol;

    A.SetInterval(0, N);
    b.SetInterval(0, N);
    sol.SetInterval(0, N);

    for(int k = 0; k < N; ++k){
        int i = k / n, j = k % n;
        double x = h * i, y = h * j;
        if(i == 0 || j == 0 || i == n - 1 || j == n - 1){ // boundary points
            A[k][k] = 1;
            b[k] = g_bond(x, y);
        } else {
            A[k][k] = 2 * (dx + dy);
            A[k][k + n] = -dx;
            A[k][k - n] = -dx;
            A[k][k - 1] = -dy;
            A[k][k + 1] = -dy;
            b[k] = f(x, y) * h * h;
        }
    }

    std::cout<< "matrix solver started\n";

    Solver S(Solver::INNER_ILU2);
    S.SetParameter("absolute_tolerance", "1e-12");
    S.SetParameter("verbosity", "2");
    S.SetParameter("relative_tolerance", "1e-7");
    S.SetMatrix(A);
    bool solved = S.Solve(b, sol);
    std::cout << "lin.it.: " << S.Iterations() << std::endl;
    if(!solved){
        std::cout << "Linear solver failure!" << std::endl;
        std::cout << "Reason: " << S.ReturnReason() << std::endl;
    }
    std::cout << "matrix solver ended\n";

    std::vector<double> full_sol(N);
    for(int i = 0; i < N; ++i){
        full_sol[i] = sol[i];
    }

    std::vector<double> exact_sol(N);
    for(int k = 0; k < N; ++k){
        int i = k / n, j = k % n;
        exact_sol[k] = u(i * h, j * h, dx, dy);
    }


    std::string res_fname = "./res.vtk";
    if(argc > 2){
        res_fname = std::string(argv[2]);
    }

    double l2_norm = 0., c_norm = 0;
    for(int i = 0; i < N; ++i){
        l2_norm += (exact_sol[i] - sol[i]) * (exact_sol[i] - sol[i]);
        c_norm = std::max(c_norm, abs(exact_sol[i] - sol[i]));
    }
    l2_norm = sqrt(l2_norm) / N;
    std::cout << "L2 norm: " << l2_norm << " C-norm: " << c_norm << std::endl;


    std::cout << "saving solution started with res in " << res_fname << " file\n";
    save_solution(n, full_sol, exact_sol, res_fname);
    std::cout << "saving solution ended\n";

    Solver::Finalize();
}