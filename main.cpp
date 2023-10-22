//
// Created by nebil on 22.10.23.
//
#include <cstdio>
#include <vector>
#include <stdexcept>
#include <iostream>

#include "mesh.h"
#include "inmost.h"

using namespace INMOST;

extern double f(double x, double y);

double u(double x, double y, double dx = 1., double dy = 1.);


void save_solution(const SquareMesh& mesh, const std::vector<double>& sol, const std::vector<double>& exact_sol, std::string vtk_fname){
    FILE* fd = fopen(vtk_fname.c_str(), "w");
    if(!fd){
        throw std::runtime_error("can't open file " + vtk_fname);
    }
    fprintf(fd, "# vtk DataFile Version 2.0\n");
    fprintf(fd, "solution by FDM\n");
    fprintf(fd, "ASCII\nDATASET UNSTRUCTURED_GRID\n");
    fprintf(fd, "POINTS %lu double\n", mesh.points.size());
    for(auto &p: mesh.points){
        fprintf(fd, "%lf %lf %lf\n", p[0], p[1], 0.0);
    }
    fputc('\n', fd);
    fprintf(fd, "CELLS %lu %lu\n", mesh.squares.size(),  mesh.squares.size() * 5);
    for(auto &it: mesh.squares){
        fprintf(fd, "4 %lu %lu %lu %lu\n", it[0], it[1], it[2], it[3]);
    }
    fputc('\n', fd);
    fprintf(fd, "CELL_TYPES %lu\n", mesh.squares.size());
    for(int i = 0; i < mesh.squares.size(); ++i){
        fprintf(fd, "9\n");
    }
    fprintf(fd, "CELL_DATA %lu\n",  mesh.squares.size());
    fprintf(fd, "SCALARS U double 1\n");
    fprintf(fd, "LOOKUP_TABLE default\n");
    for(int i = 0; i < mesh.squares.size(); ++i){
        fprintf(fd, "%lf\n", sol[i]);
    }
    fputc('\n', fd);

    if(!exact_sol.empty()) {
        fprintf(fd, "SCALARS U_true double 1\n");
        fprintf(fd, "LOOKUP_TABLE default\n");
        for (int i = 0; i < mesh.squares.size(); ++i) {
            fprintf(fd, "%lf\n", exact_sol[i]);
        }
        fputc('\n', fd);

        fprintf(fd, "SCALARS rtol double 1\n");
        fprintf(fd, "LOOKUP_TABLE default\n");
        for (int i = 0; i < mesh.squares.size(); ++i) {
            fprintf(fd, "%lf\n", abs(exact_sol[i] - sol[i]) / abs(exact_sol[i]));
        }
        fputc('\n', fd);
    }
}


int main(int argc, char**argv){
    Solver::Initialize(&argc, &argv);
    std::cout << "parsing " << std::string(argv[1]) << " file\n";
    SquareMesh mesh(argv[1]);
    std::cout << "h: " << mesh.h() <<" n: "<< mesh.size().first <<" m: "<< mesh.size().second << '\n';
    std::cout << "cells: " << mesh.squares.size() << std::endl;



    size_t N = mesh.squares.size(), n = mesh.size().first, m = mesh.size().second;
    double h = mesh._h;

    double dx = 4.0, dy = 1.0;

    Sparse::Matrix A;
    Sparse::Vector b;
    Sparse::Vector sol;

    A.SetInterval(0, N);
    b.SetInterval(0, N);
    sol.SetInterval(0, N);

    for(int k = 0; k < N; ++k){
        A[k][k] = 2 * (dx + dy);
        auto temp = mesh.index_2d(k);
        size_t i = temp.first, j = temp.second;
        if(i < N - 1){
            A[k][k + n] = -dx;
        }
        if(i > 0){
            A[k][k - n] = -dx;
        }
        if(j > 0){
            A[k][k - 1] = -dy;
        }
        if(j < N - 1){
            A[k][k + 1] = -dy;
        }
        point p = mesh.centers[k];
        b[k] = f(p[0], p[1]) * h * h;
    }

    std::cout<< "matrix solver started\n";

    Solver S(Solver::INNER_ILU2);
    S.SetParameter("absolute_tolerance", "1e-10");
    S.SetParameter("verbosity", "2");
    S.SetParameter("relative_tolerance", "1e-6");
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
    for(int i = 0; i < N; ++i){
        point p = mesh.centers[i];
        exact_sol[i] = u(p[0], p[1], dx, dy);
    }


    std::string res_fname = "./res.vtk";
    if(argc > 2){
        res_fname = std::string(argv[2]);
    }

    std::cout << "saving solution started with res in " << res_fname << " file\n";
    save_solution(mesh, full_sol, exact_sol, res_fname);
    std::cout << "saving solution ended\n";
    Solver::Finalize();
}