#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>

#include "progressbar.h"

using std::vector;
using namespace std::chrono;


/**
 * Apply rotation matrix to A (rotating 'i' and 'j') rows
 * @param A original matrix
 * @param i first row
 * @param j second row
 * @param c cos of rotation
 * @param s sin of rotation
 * @param start first != 0 element of rows.
 */
void rotationRow(vector<vector<double>>& A, int i, int j, double c, double s, int start = 0){
    size_t n = A.size();
    for(int k = start; k < n; ++k){
        double aik = c * A[i][k] - s * A[j][k];
        double ajk = s * A[i][k] + c * A[j][k];
        A[i][k] = aik, A[j][k] = ajk;
    }
}

/**
 * Apply rotation matrix to A (rotation 'i' and 'j') cols
 * @param A original matrix
 * @param i first col
 * @param j second col
 * @param c cos of rotation
 * @param s sin of rotation
 * @param start first != 0 element of cols
 */
void rotationCol(vector<vector<double>>& A, int i, int j, double c, double s, int start = 0){
    size_t n = A.size();
    for(int k = start; k < n; ++k){
        double aki = c * A[k][i] - s * A[k][j];
        double akj = s * A[k][i] + c * A[k][j];
        A[k][i] = aki, A[k][j] = akj;
    }
}

void print_matrix(const vector<vector<double>>& A){
    size_t n = A.size();
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            printf("%4lf ", A[i][j]);
        }
        putchar('\n');
    }
    putchar('\n');
}

/**
 * Tridiagonalize of original matrix A.
 * @param A original matrix.
 */
void Tridiagonaliation(vector<vector<double>>&A){
    size_t n = A.size();
    progressbar bar(n);
    for(int i = 0; i < n; ++i){
        bar.update();
        for(int j = i + 2; j < n; ++j){
            double mod = sqrt(A[i+1][i] * A[i+1][i] + A[j][i] * A[j][i]);
            double s = -A[j][i] / mod, c = A[i+1][i] / mod;
            rotationRow(A,i+1, j, c, s,i);
            rotationCol(A,i+1, j, c, s,i);
        }
    }
}


int main() {
    int n;
    scanf("%d", &n);

    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<double> distr(-100, 100);

    std::vector<std::vector<double>> A(n, std::vector<double>(n));
    for(int i = 0; i < n; ++i){
        for(int j = i; j < n; ++j){
            A[i][j] = distr(eng);
            A[j][i] = A[i][j];
        }
    }
    //print_matrix(A);
    auto start = high_resolution_clock::now();
    Tridiagonaliation(A);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout <<std::endl<< "duration: "<< duration.count() << std::endl;
    //print_matrix(A);
}
