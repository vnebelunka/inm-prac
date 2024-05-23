//
// Created by nebil on 22.10.23.
//
#include <cmath>

double f(double x, double y){
    return sin(10 * x) * sin(10 * y);
}

double u(double x, double y, double dx = 1., double dy = 1.){
    return sin(10 * x) * sin(10 * y) / ((dx + dy) * 10 * 10);
}

double g_bond(double x, double y){
    return sin(10 * x) * sin(10 * y) / ((2) * 10 * 10);
}