/**
 * Problem：243 Runge-Kutta方法
 * Time：2024.07.01
 */

#include <stdio.h>
#include <stdlib.h>

double *RK4(double (*func)(double, double), double a, double b, double h, double y0); // 经典四阶四级Runge-Kutta方法

// y'=f(x,y)
double f(double x, double y)
{
    return y - 2 * x / y;
}

int main()
{
    double a = 0, b = 1, h = 0.1, y0 = 1, *y;
    int n = 11;
    y = RK4(f, a, b, h, y0);
    for (int i = 0; i < n; ++i)
    {
        printf("y[%d] = %lf\n", i, y[i]);
    }
    free(y);

    return 0;
}

double *RK4(double (*func)(double, double), double a, double b, double h, double y0)
{
    int n = (int)((b - a) / h) + 1;
    double *x = (double *)malloc(n * sizeof(double));
    double *y = (double *)malloc(n * sizeof(double));

    // 赋初值
    for (int i = 0; i < n; ++i)
        x[i] = a + i * h;
    y[0] = y0;

    // 迭代格式
    double K0, K1, K2, K3;
    for (int i = 1; i < n; i++)
    {
        K0 = func(x[i - 1], y[i - 1]);
        K1 = func(x[i - 1] + h / 2, y[i - 1] + K0 * h / 2);
        K2 = func(x[i - 1] + h / 2, y[i - 1] + K1 * h / 2);
        K3 = func(x[i], y[i - 1] + K2 * h);
        y[i] = y[i - 1] + h * (K0 + 2 * K1 + 2 * K2 + K3) / 6;
    }

    return y;
}