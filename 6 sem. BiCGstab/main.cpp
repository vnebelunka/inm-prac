#include "Matrix.hpp"
#include <random>

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

template <typename T>
vector<T> BiCGstab(Matrix<T> A, vector<T> b, double eps){
    int n = b.size();
    vector<T> x(n), r = b - matvec(A, x), rconj = r, p = r;
    int iter = 0;
    for(;iter < 1e6; ++iter){
        vector<T> Ap = matvec(A, p);
        T a = dot(r, rconj) / dot(Ap, rconj);
        vector<T> s = r - a * Ap;
        if(fabs(norm(s)) < eps){
            x = x + a * p;
            break;
        }
        vector<T> As = matvec(A, s);
        T w = dot(As, s) / dot(As, As);
        x = x + a * p + w * s;
        T bk = 1. / dot(r, rconj);
        // orthogonality check
        if(iter % 1000 != 1) {
            r = s - w * As;
        } else {
            r = b - matvec(A, x);
        }
        if(fabs(norm(r)) < eps){
            break;
        }
        bk *= dot(r, rconj);
        p = r + bk * p - w * Ap;
        // restart
        if(fabs(dot(r, rconj)) < 1e-8){
            rconj = r;
            p = r;
        }
    }
    std::cout << "iters:" << iter << std::endl;
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
    std::cout << "x = \n";
    std::cout <<"r = " << norm(b - matvec(A, x)) << std::endl;
}