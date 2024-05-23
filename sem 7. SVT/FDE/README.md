# Постановка задачи
Решается двумерная задача Дирихле для двумерного стационарного оператора диффузии
в единичном квадрате на регулярной квадратной сетке методом элементов.

Подробнее - в report.pdf

# Как запускать
```console
./FDM [mesh.vtk]
```
В файле _function.hpp_ находятся функции, задающие задачу, а также тензор диффузии.

```C++
const double dx = 1.0;
const double dy = 1.0;
const double dxy = 0.0;

double source(double x, double y); // правая часть

double conc_an(double x, double y, double dx = 1., double dy = 1.); // аналитическое решение

```