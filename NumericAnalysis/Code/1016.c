/**
 * Problem:1016 割线法求根
 * Time:2024.03.22
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double secant(double x1, double x2); // 割线法求根函数
double f(double x);                  // 待求根的函数 F(x)

int main()
{
    double x1, x2; // 区间左右端点

    scanf("%lf%lf", &x1, &x2);
    printf("%lf", secant(x1, x2));
}

double f(double x)
{
    return cos(x) - tan(x);
}

double secant(double x1, double x2)
{
    while (fabs(x1 - x2) > 1e-8)
    { // 循环结束要求：区间长小于1e-8
        double tmp = x1 - f(x1) * (x2 - x1) / (f(x2) - f(x1));
        x1 = x2;
        x2 = tmp;
    }

    return x1;
}