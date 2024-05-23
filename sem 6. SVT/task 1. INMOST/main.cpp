#include <iostream>
#include "inmost.h"

using namespace INMOST;

int main(int argc, char**argv){
    Solver::Initialize(&argc, &argv);
    Sparse::Matrix A;
    double tempTimer, solvingTimer;
    bool success;

    char *matrix_fname = argv[1], *rhs_fname = argv[2], *drop_tol = argv[3];

    A.Load(matrix_fname);
    Sparse::Vector rhs;
    rhs.Load(rhs_fname);

    Solver S("inner_ilu2");

    S.SetParameter("verbosity", "2");
    S.SetParameter("drop_tolerance", drop_tol);

    tempTimer = Timer();
    S.SetMatrix(A);
    std::cout << "preconditioner time: " << Timer() - tempTimer << std::endl;
    Sparse::Vector sol(rhs);
    tempTimer = Timer();
    success = S.Solve(rhs, sol);
    solvingTimer = Timer() - tempTimer;
    std::cout << "solving time: " << solvingTimer << std::endl;
    std::cout << "iterations: " << S.Iterations() << std::endl;
    Solver::Finalize();
    return 0;
}