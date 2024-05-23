//
// Created by nebil on 22.10.23.
//
#include <cmath>

double dx = 7.75, dy = 3.25, dxy = 3.89711432;

//double dx = 10, dy = 1, dxy = 0;

const double a = 10;

double u(double x, double y){
    return sin(a*x) * sin(a*y);
}

double f(double x, double y){
    return -a*a * (2.*dxy * cos(a*x)*cos(a*y) - (dx+dy) * sin(a*x)*sin(a*y));
}