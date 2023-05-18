#pragma once

#include <vector>

using std::vector;

template<typename matrix>
vector<double> BiCGstab(const matrix &A_chunk, vector<double>& b_chunk, double rel_tol);

#include "BiCGstab.tpp"