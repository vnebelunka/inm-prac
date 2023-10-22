//
// Created by nebil on 22.10.23.
//
#include <cmath>

double f(double x, double y){
    return sin(M_PI * x) * sin(M_PI * y);
}

double u(double x, double y, double dx = 1., double dy = 1.){
    return sin(M_PI * x) * sin(M_PI * y) / ((dx + dy) * M_PI * M_PI);
}

