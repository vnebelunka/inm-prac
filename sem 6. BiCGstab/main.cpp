#include "Matrix.hpp"
#include <random>
#include <complex>

using std::vector;
using std::complex;

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
                A(i, i) += fabs(A(i, j));
            }
        }
    }
}

template <typename DataType>
vector<DataType> BiCGstab(Matrix<DataType> A, vector<DataType> b, double eps){
    int n = b.size();
    vector<DataType> x(n), r_k = b - matvec(A, x), rconj0 = r_k, p = r_k;
    int iter = 0;
    for(;iter < 1e4; ++iter){
        vector<DataType> Ap_k = matvec(A, p);
        DataType alpha_k = dot(r_k, rconj0) / dot(Ap_k, rconj0);
        vector<DataType> s_k = r_k - alpha_k * Ap_k;
        if(fabs(norm(s_k)) < eps){
            x = x + alpha_k * p;
            break;
        }
        vector<DataType> As_k = matvec(A, s_k);
        DataType w_k = dot(As_k, s_k) / dot(As_k, As_k);
        x = x + alpha_k * p + w_k * s_k;
        DataType r_krconj_0 = dot(r_k, rconj0);
        // orthogonality check
        if(iter % 1000 != 1) {
            r_k = s_k - w_k * As_k;
        } else {
            r_k = b - matvec(A, x);
            std::cout << dot(r_k, r_k) << std::endl;
        }
        if(fabs(norm(r_k)) < eps){
            break;
        }
        DataType b_k = dot(r_k, rconj0) / r_krconj_0;
        p = r_k + b_k * p - w_k * Ap_k;
        // restart
        if(fabs(dot(r_k, rconj0)) < 1e-8){
            rconj0 = r_k;
            p = r_k;
        }
    }
    std::cout << "iters: " << iter << std::endl;
    return x;
}


int main(){
    int n;
    std::cin >> n;
    Matrix<double> A(n, n);
    vector<double> b(n);
    gen_matrix(A);
    if(n < 10) {
        print_matrix(A);
    }
    for(int i = 0; i < n; ++i){
        b[i] = sin(i);
    }

    vector<double> x = BiCGstab(A, b, 1e-8);
    std::cout <<"r = " << norm(b - matvec(A, x)) << std::endl;
}
