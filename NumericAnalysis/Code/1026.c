/**
 * Problem:1026 Runge-Kutta方法
 * Time：2024.04.28
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define E 2.718281828459

double f(double x, double y);
void partition(int n, double left, double right, double x[]); // 分划区间函数
void RungeKutta2(int n, double x[], double y[]);              // 二级 RungeKutta 方法
void RungeKutta4(int n, double x[], double y[]);              // 四级 RungeKutta 方法

int main()
{

    int n = 65; // 分点个数（ = 小区间个数 + 1）
    double left = 0, right;

    double *x = (double *)malloc(sizeof(double) * n); // 利用 malloc 函数为数组 x,y 申请空间
    double *y = (double *)malloc(sizeof(double) * n);

    scanf("%lf", &right);

    partition(n, left, right, x);
    y[0] = 0;
    RungeKutta2(n, x, y);
    RungeKutta4(n, x, y);

    // 输出精确解
    printf("%lf", right * pow(E, -right));

    // 释放数组 x,y
    free(x);
    free(y);

    return 0;
}

double f(double x, double y)
{
    return pow(E, -x) - y;
}

void partition(int n, double left, double right, double x[])
{
    double h = (right - left) / (n - 1); // 区间长
    for (int i = 0; i <= n - 1; i++)
    {
        // left = x[0] < x[1] < ... < x[n - 1] = right
        x[i] = left + h * i;
    }
    return;
}

void RungeKutta2(int n, double x[], double y[])
{
    // 二级 Runge-Kutta 方法
    double k1, k2, h = x[1] - x[0];
    for (int i = 0; i <= n - 2; i++)
    {
        k1 = f(x[i], y[i]);
        k2 = f(x[i] + h / 2, y[i] + h * k1 / 2);
        y[i + 1] = y[i] + h * k2;
    }
    printf("%lf\n", y[n - 1]);
}

void RungeKutta4(int n, double x[], double y[])
{
    // 四级 Runge-Kutta 方法
    double k1, k2, k3, k4, h = x[1] - x[0];
    for (int i = 0; i <= n - 2; i++)
    {
        k1 = f(x[i], y[i]);
        k2 = f(x[i] + h / 2, y[i] + h * k1 / 2);
        k3 = f(x[i] + h / 2, y[i] + h * k2 / 2);
        k4 = f(x[i] + h, y[i] + h * k3);
        y[i + 1] = y[i] + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    }
    printf("%lf\n", y[n - 1]);
}
