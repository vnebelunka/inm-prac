# Постановка задачи
Решается двумерная задача Дирихле для двумерного нестационарного оператора диффузии
в единичном квадрате на регулярной квадратной сетке с дискретизацией методом конечных разностей (явным и неявным методом).

Подробнее - в report.pdf

# Как запускать
```console
./FDM n 
```
*  n - колчичество точек в разбиении стороны квадрата

В файле _function.hpp находятся функции, задающие задачу
```C++
enum time_iteration_type
{
	EXPLICIT = 1,
	IMPLICIT = 2
};

const time_iteration_type Problem_type = EXPLICIT; // метод решения по времени
const double dx = 7.75, dy = 3.25, dxy = 3.89711432; // параметры тензора диффузии
const double dt = 3e-4; // шаг по времени
const int T = 100; // кол-во шагов по времени


double u0(double x, double y){ // начальное условие
    if(x > 0.35 && x < 0.65 && y > 0.35 && y < 0.65){
        return 1;
    }
    return 0;
}

double phi(double x, double y, double t){ // граничное условие
    return 0;
}

double f(double x, double y, double t){ // функция - источник
    return 0;
}
```