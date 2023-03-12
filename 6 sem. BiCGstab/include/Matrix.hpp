#pragma once

#include <vector>
#include <iostream>
template<class T>
class Matrix {
    std::vector<T> data_;
    size_t n,  m;
public:
    Matrix(size_t n, size_t m = 1): n(n), m(m){
        data_.resize(n * m);
    }

    T& operator()(int i, int j){
        return data_[n * i + j];
    }

    T operator()(int i, int j) const{
        return  data_[n * i + j];
    }

    std::pair<size_t, size_t> size() const{
        return {n, m};
    }
};

void print_matrix(Matrix<double> const& A);

template<typename T>
T dot(const std::vector<T>& a,const std::vector<T>& b);

template<typename T>
std::vector<T> matvec(const Matrix<T>& A, const std::vector<T>& x);

template <typename T>
std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &b);

template <typename T>
std::vector<T> operator*(const std::vector<T> &a, T b);

template <typename T>
std::vector<T> operator*(T b, const std::vector<T> &a);

template <typename T>
std::vector<T> operator-(const std::vector<T> &a, const std::vector<T> &b);

template <typename T>
double norm(const std::vector<T>& a);

#include "Matrix.tpp"