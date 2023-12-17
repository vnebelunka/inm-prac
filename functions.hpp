//
// Created by nebil on 22.10.23.
//
#include <cmath>

const double dx = 5.0;
const double dy = 1.0;
const double dxy = 0.0;
const double a = 4;

double C(double x, double y) // analytical solution
{
	//return x;
	return sin(a*x) * sin(a*y);
}

double source(double x, double y) 
{
	//return 0;
	return -a*a * (2.*dxy * cos(a*x)*cos(a*y) - (dx+dy) * sin(a*x)*sin(a*y));
}