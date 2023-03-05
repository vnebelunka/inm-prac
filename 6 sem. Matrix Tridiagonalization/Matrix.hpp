#include <vector>
#include <iostream>

#ifndef CUCUMBER_MATRIX_HPP
#define CUCUMBER_MATRIX_HPP
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


#endif //CUCUMBER_MATRIX_HPP
