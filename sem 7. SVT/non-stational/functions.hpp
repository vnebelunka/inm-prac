//
// Created by nebil on 22.10.23.
//
#include <cmath>

enum time_iteration_type
{
	EXPLICIT = 1,
	IMPLICIT = 2
};

const time_iteration_type Problem_type = EXPLICIT;
const double dx = 7.75, dy = 3.25, dxy = 3.89711432;
const double dt = 3e-4;
const int T = 100;


double u0(double x, double y){
    if(x > 0.35 && x < 0.65 && y > 0.35 && y < 0.65){
        return 1;
    }
    return 0;
}

double phi(double x, double y, double t){
    return 0;
}

double f(double x, double y, double t){
    return 0;
}