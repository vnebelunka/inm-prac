#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include "Matrix.hpp"


using namespace std::chrono;

extern "C" {
    extern void dsytd2_(char* UPLO, int* N, double*A, int* LDA, double* D, double *E, double *tau, int *info);
}

int main() {
    int n;
    scanf("%d", &n);

    //std::random_device rd;
    std::default_random_engine eng;
    std::uniform_real_distribution<double> distr(-100, 100);

    Matrix<double> A(n, n);
    for(int i = 0; i < n; ++i){
        for(int j = i; j < n; ++j){
            A(i, j) = distr(eng);
            A(j, i) = A(i, j);
        }
    }
    double *D = (double*) calloc(n, sizeof(double));
    double *E = (double*) calloc(n-1, sizeof(double));
    double *tau = (double*) calloc(n-1, sizeof(double));
    int info = 0;
    char UPLO = 'U';
    int LDA = n;
    auto start = high_resolution_clock::now();
    dsytd2_(&UPLO, &n, A.data(), &LDA, D, E, tau, &info);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout <<std::endl<< "duration: "<< duration.count() << std::endl;
    free(D), free(E), free(tau);
}
