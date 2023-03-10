#include <iostream>
#include "inmost.h"

using namespace INMOST;

int main(){
    Sparse::Matrix A;
    A.Load("/home/vnebelunka/Загрузки/2/A_tet1.mtx");
    Sparse::Vector rhs;
    rhs.Load("/home/vnebelunka/Загрузки/2/rhs_tet1.mtx");

    Solver S("inner_ilu2");
    S.SetParameter("verbosity", "3");
    S.SetMatrix(A);
    Sparse::Vector sol(rhs);
    std ::cout  << S.Solve(sol, rhs) << std::endl;
    return 0;
}