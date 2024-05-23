#pragma once

#include <complex.h>
#include "Matrix.hpp"
#include <cblas.h>
#include <omp.h>

template<typename T>
T dot(const std::vector<T>& a,const std::vector<T>& b){
    T ans = 0;
#pragma omp parallel for shared(ans, a, b) default(none)
    for(int i = 0; i < a.size(); ++i){
        ans += a[i] * b[i];
    }
    return ans;
}

using std::complex;



template<typename T>
complex<T> dot(const std::vector<complex<T>>& a, const std::vector<complex<T>>& b){
    complex<T> ans = 0;
#pragma omp parallel for shared(ans, a, b) default(none)
    for(int i = 0; i < a.size(); ++i){
        ans += a[i] * conj(b[i]);
    }
    return ans;
}

template<typename T>
std::vector<T> matvec(const Matrix<T> & A, const std::vector<T>& x){
    size_t n = A.size().first, m = A.size().second;
    std::vector<T> ans(m);
#pragma omp parallel for shared(ans, x, n, m, A) default(none) collapse(2)
    for(int i = 0; i < m; ++i){
        for(int j = 0; j < n; ++j){
            ans[i] += A(i, j) * x[j];
        }
    }
    return ans;
}

template <typename T>
std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &b){
    std::vector<T> ans(a.size());
#pragma omp parallel for shared(ans, a, b) default(none)
    for(int i = 0; i < a.size(); ++i){
        ans[i] = a[i] + b[i];
    }
    return ans;
}

template <typename T>
std::vector<T> operator-(const std::vector<T> &a, const std::vector<T> &b){
    std::vector<T> ans(a.size());
#pragma omp parallel for shared(ans, a, b) default(none)
    for(int i = 0; i < a.size(); ++i){
        ans[i] = a[i] - b[i];
    }
    return ans;
}

template <typename T>
std::vector<T> operator*(const std::vector<T> &a, T b){
    std::vector<T> ans(a.size());
#pragma omp parallel for shared(ans, a, b) default(none)
    for(int i = 0; i < a.size(); ++i){
        ans[i] = a[i] * b;
    }
    return ans;
}

template <typename T>
std::vector<T> operator*(T b, const std::vector<T> &a){
    return a * b;
}

template <typename T>
double norm(const std::vector<T>& a){
    double ans = 0;
#pragma omp parallel for shared(ans, a) default(none)
    for(int i = 0; i < a.size(); ++i){
        ans += a[i] * a[i];
    }
    return sqrt(ans);
}

template <typename T>
double norm(const std::vector<std::complex<T>>& a){
    double ans = 0;
#pragma omp parallel for shared(ans, a) default(none)
    for(int i = 0; i < a.size(); ++i){
        ans += a[i] * conj(a[i]);
    }
    return sqrt(ans);
}
