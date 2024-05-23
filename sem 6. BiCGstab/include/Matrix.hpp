#pragma once

#include <vector>
#include <iostream>
#include <mpi.h>


template<class T>
class Matrix {
    std::vector<T> data_;
    size_t n,  m;
    int proc_id, n_proc;
    MPI_Comm comm;
public:

    Matrix(int proc_id=0, int n_proc=1, MPI_Comm comm=MPI_COMM_WORLD): proc_id(proc_id), n_proc(n_proc), comm(comm){n = 0, m = 0;}

    Matrix(size_t n, size_t m = 1): n(n), m(m){
        data_.resize(n * m);
    }

    void resize(size_t n, size_t m){
        data_.resize(n * m);
        this->n = n, this->m = m;
    }

    T& operator()(int i, int j){
        return data_[m * i + j];
    }

    T operator()(int i, int j) const{
        return  data_[m * i + j];
    }

    std::pair<size_t, size_t> size() const{
        return {n, m};
    }

    T* data(){
        return data_.data();
    }

    std::vector<T> matvec(const std::vector<T> &x_chunk) const;

    int get_n_proc() const{return n_proc;}
    int get_proc_id() const{return proc_id;}
    MPI_Comm get_comm() const {return comm;}
};

void print_matrix(Matrix<double> const& A);
void print_matrix(Matrix<double> const& A, int proc_id);

template<typename T>
T dot(const std::vector<T>& a,const std::vector<T>& b);

double dot( std::vector<double> &a_chunk, std::vector<double> & b_chunk, int proc_id, int n_proc, MPI_Comm comm);

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

double norm(std::vector<double> &a_chunk, int proc_id, int n_proc, MPI_Comm comm);

#include "Matrix.tpp"