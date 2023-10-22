//
// Created by nebil on 22.10.23.
//

#ifndef FDM_MESH_H
#define FDM_MESH_H
#include <vector>
#include <array>

typedef std::array<double, 2> point;

struct SquareMesh{
    using point = std::array<double, 2>;

    std::vector<point> points;
    std::vector<std::array<size_t, 4>> squares;
    std::vector<point> centers;
    double _h;
    int _n, _m;

    SquareMesh(const char* vtk_fname);
    double h() const {return _h;}
    std::pair<size_t, size_t> size() const {return {_n, _m};}
    std::pair<size_t, size_t> index_2d(size_t idx) const {return {idx / _m, idx % _m};}
    size_t index_1d(size_t i, size_t j) const{return i * _m + j;}
    size_t index_1d(std::pair<size_t, size_t> idx) const{return idx.first * _m + idx.second;}
    point center(size_t i, size_t j) const{return centers[index_1d(i, j)];}
    point center(std::pair<size_t, size_t> idx) const{return centers[index_1d(idx)];}
};

#endif //FDM_MESH_H
