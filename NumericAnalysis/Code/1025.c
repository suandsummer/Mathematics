/**
 * Problem:1025 简单数值积分法
 * Time:2024.04.28
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f(double x, double y);                                  // f(x,y)
void forwardEuler(int n, double x[], double y[]);              // 向前Euler方法
void backwardEuler(int n, double x[], double y[], double eps); // 向后Euler方法
void trapezoid(int n, double x[], double y[], double eps);     // 梯形方法
void improvedEuler(int n, double x[], double y[]);             // 改进Euler方法
void partition(int n, double left, double right, double x[]);  // 区间划分

int main()
{

    int n = 5; // 分划节点个数

    // 利用 malloc 函数为数组 x[],y[]开辟内存空间以记录分划端点信息
    double *x = (double *)malloc(sizeof(double) * n);
    double *y = (double *)malloc(sizeof(double) * n);

    y[0] = 1; // 设定初值 y(0) = 1

    double left = 0, right; // 左右区间端点
    double eps = 1.0e-10;   // 迭代容许误差

    scanf("%lf", &right);

    // 对区间做等距分划
    partition(n, left, right, x);

    // 调用函数，计算微分方程的解并输出
    forwardEuler(n, x, y);
    backwardEuler(n, x, y, eps);
    trapezoid(n, x, y, eps);
    improvedEuler(n, x, y);

    // 输出精确解
    printf("%lf\n", sqrt(1 + 2 * right));

    // 利用 free 函数释放数组 x[] , y[]
    free(x);
    free(y);

    return 0;
}

double f(double x, double y)
{

    return y - 2 * x / y; // 在该实例中， f(x,y) 取为 y - 2x / y
}

void partition(int n, double left, double right, double x[])
{

    double h = (right - left) / (n - 1); // 区间步长

    // 记录分划端点信息到数组 x[] 中
    for (int i = 0; i <= n - 1; i++)
        x[i] = left + h * i;

    return;
}

void forwardEuler(int n, double x[], double y[])
{
    // 向前 Euler 方法
    double h = x[1] - x[0];
    for (int i = 1; i < n; ++i)
        y[i] = y[i - 1] + h * f(x[i - 1], y[i - 1]);
    // 输出向前 Euler 方法计算得到的答案
    printf("%lf\n", y[n - 1]);

    return;
}

void backwardEuler(int n, double x[], double y[], double eps)
{
    // 向后 Euler 方法
    double f1, f2, last, err;
    double h = x[1] - x[0];
    for (int i = 1; i < n; ++i)
    {
        f1 = f(x[i - 1], y[i - 1]);
        y[i] = y[i - 1] + h * f1;
        err = 1 + eps;
        while (err >= eps)
        {
            last = y[i];
            f2 = f(x[i], y[i]);
            y[i] = y[i - 1] + h * f2;
            err = fabs(y[i] - last);
        }
    }
    // 输出向后 Euler 方法计算得到的答案
    printf("%lf\n", y[n - 1]);

    return;
}

void trapezoid(int n, double x[], double y[], double eps)
{
    // 梯形方法
    double f1, f2, tmp, err = 1 + eps;
    double h = x[1] - x[0];
    for (int i = 1; i < n; ++i)
    {
        f1 = f(x[i - 1], y[i - 1]);
        y[i] = y[i - 1] + h * f1;
        err = 1 + eps;
        while (err > eps)
        {
            f2 = f(x[i], y[i]);
            tmp = y[i - 1] + h * (f1 + f2) / 2;
            err = fabs(y[i] - tmp);
            y[i] = tmp;
        }
    }
    // 输出梯形方法计算得到的答案
    printf("%lf\n", y[n - 1]);

    return;
}

void improvedEuler(int n, double x[], double y[])
{
    // 改进 Euler 方法
    double h = x[1] - x[0];
    for (int i = 1; i < n; ++i)
    {
        double f1 = f(x[i - 1], y[i - 1]);
        y[i] = y[i - 1] + h * f1;
        double f2 = f(x[i], y[i]);
        y[i] = y[i - 1] + h * (f1 + f2) / 2;
    }
    // 输出改进 Euler 方法计算得到的答案
    printf("%lf\n", y[n - 1]);

    return;
}
