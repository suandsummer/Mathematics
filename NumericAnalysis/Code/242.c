/**
 * Problem：242 向前欧拉法
 * Time：2024.07.01
 */

#include <stdio.h>
#include <stdlib.h>

double *ForEuler(double (*func)(double, double), double a, double b, double h, double y0); // 向前欧拉法
double f(double x, double y);                                                              // y'=f(x,y)

int main()
{
    double a = 0, b = 1, h = 0.1, y0 = 1, *y;
    int n = 11;
    y = ForEuler(f, a, b, h, y0);
    for (int i = 0; i < n; ++i)
    {
        printf("y[%d] = %lf\n", i, y[i]);
    }
    free(y);

    return 0;
}

double f(double x, double y)
{
    return y - 2 * x / y;
}

double *ForEuler(double (*func)(double, double), double a, double b, double h, double y0)
{
    int n = (int)((b - a) / h) + 1;
    double *y = (double *)malloc(n * sizeof(double));
    y[0] = y0;
    for (int i = 1; i < n; i++)
        y[i] = y[i - 1] + h * func(a + (i - 1) * h, y[i - 1]);
    return y;
}