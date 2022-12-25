#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>

#include <cblas.h>
#include "Matrix.hpp"


using std::vector;
using namespace std::chrono;

extern "C" {
    extern void dsytd2(char UPLO, int N, double*A,int LDA, double*D, double *E, double*tau, int *info);
}



/**
 * Apply rotation matrix to A (rotating 'i' and 'j') rows
 * @param A original matrix
 * @param i first row
 * @param j second row
 * @param c cos of rotation
 * @param s sin of rotation
 * @param start first != 0 element of rows.
 */
void rotationRow(Matrix<double> &A, int i, int j, double c, double s, int start = 0){
    size_t n = A.size().first;
    for(int k = start; k < n; ++k){
        double aik = c * A(i, k) - s * A(j, k);
        double ajk = s * A(i, k) + c * A(j, k);
        A(i, k) = aik, A(j, k) = ajk;
    }
}


/**
 * Checks if A == Q^T D Q
 * @param A
 * @param D
 * @param Q
 */
void check(Matrix<double> const &A, Matrix<double> const &D, vector<std::pair<double, double>> &rots){
    size_t n = A.size().first;
    Matrix<double> C(n, n), F(n, n);
    Matrix<double> Q(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Q(i, j) = i == j;
        }
    }
    int pos = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 2; j < n; ++j) {
            double c = rots[pos].first, s = rots[pos++].second;
            rotationRow(Q, i + 1, j, c, s, 0);
        }
    }


    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            for(int k = 0; k < n; ++k){
                C(i, j) += Q(k, i) * D(k, j);
            }
        }
    }
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            for(int k = 0; k < n; ++k){
                F(i, j) += C(i, k) * Q(k, j);
            }
        }
    }
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            if(fabs(A(i, j) - F(i, j)) >= 1e-10){
                std::cout << "Matrix factorization is incorrect!\n";
                exit(1);
            }
        }
    }
    std::cout << "Matrix factorization is correct!\n";
}

int main() {
    int n;
    scanf("%d", &n);

    //std::random_device rd;
    std::default_random_engine eng;
    std::uniform_real_distribution<double> distr(-100, 100);

    Matrix<double> A(n, n);
    vector<std::pair<double, double>> rots;
    for(int i = 0; i < n; ++i){
        for(int j = i; j < n; ++j){
            A(i, j) = distr(eng);
            A(j, i) = A(i, j);
        }
    }
    auto start = high_resolution_clock::now();
    //Tridiagonaliation(A, D, rots);
    double *D = (double*) calloc(n, sizeof(double));
    double *E = (double*) calloc(n-1, sizeof(double));
    double *tau = (double*) calloc(n-1, sizeof(double));
    int info[0];
    dsytd2('U', n, A.data(), 1, D, E, tau, info);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout <<std::endl<< "duration: "<< duration.count() << std::endl;
    free(D), free(E), free(tau);
}
