//
// Created by nebil on 22.10.23.
//
#include <cstdio>
#include <stdexcept>

#include "mesh.h"

SquareMesh::SquareMesh(const char *vtk_fname) {
    FILE* fd = fopen(vtk_fname, "r");
    if(!fd){
        throw std::runtime_error("can't open file " + std::string(vtk_fname));
    }
    // skipping header
    char buf[100];
    fgets(buf, 100, fd);
    fgets(buf, 100, fd);
    fgets(buf, 100, fd);
    fgets(buf, 100, fd);

    // points reader
    double x_min = 1e9, y_min = 1e9;
    double x_max = -1e9, y_max = -1e9;

    size_t num_points = 0;
    fscanf(fd, "POINTS %d", &num_points);
    fgets(buf, 100, fd);
    points.reserve(num_points);
    for (int i = 0; i < num_points; ++i) {
        point p;
        double z;
        fscanf(fd, "%lf %lf %lf", &p[0], &p[1], &z);
        if(std::abs(z) > 1e-8){
            throw std::runtime_error("mesh is not in OXY plane");
        }
        x_min = std::min(x_min, p[0]);
        x_max = std::max(x_max, p[0]);
        y_min = std::min(y_min, p[1]);
        y_max = std::max(y_max, p[1]);
        points.emplace_back(p);
    }

    // squares reader
    int n_cells;
    fscanf(fd, "\nCELLS %d", &n_cells);
    fgets(buf, 100, fd);
    squares.reserve(n_cells);
    for(int i = 0; i < n_cells; ++i){
        int type;
        fscanf(fd, "%d", &type);
        if(type == 4) {
            std::array<size_t, 4> square;
            fscanf(fd, "%lu %lu %lu %lu", &square[0], &square[1], &square[2], &square[3]);
            squares.emplace_back(square);

            point center;
            center[0] = (points[square[0]][0] + points[square[1]][0] + points[square[2]][0] + points[square[3]][0]) / 4;
            center[1] = (points[square[0]][1] + points[square[1]][1] + points[square[2]][1] + points[square[3]][1]) / 4;
            centers.push_back(center);
        } else {
            fgets(buf, 100, fd);
        }
    }
    _h = std::max(std::abs(points[squares[0][0]][0] - points[squares[0][1]][0]), std::abs(points[squares[0][0]][1] - points[squares[0][1]][1]));
    _m = ((x_max - x_min) / _h);
    _n = ((y_max - y_min) / _h);
}