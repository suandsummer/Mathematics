/**
 * Problem:1019 三次Hermite插值
 * Time:2024.03.27
 */

#include <stdio.h>
#include <math.h>

double Hermite3(double x[], double y[], double m[], double xs);

int main()
{
    double x[2], y[2], m[2];
    double xs;

    // 输入插值节点
    for (int i = 0; i < 2; i++)
    {
        scanf("%lf %lf %lf ", &x[i], &y[i], &m[i]);
    }
    scanf("%lf", &xs);

    printf("%lf", Hermite3(x, y, m, xs));

    return 0;
}

double Hermite3(double x[], double y[], double m[], double xs)
{
    double a0 = (xs - x[1]) / (x[0] - x[1]); // 插值基函数
    double a1 = (xs - x[0]) / (x[1] - x[0]);

    double b0 = 1 / (x[0] - x[1]);
    double b1 = 1 / (x[1] - x[0]);

    double res = 0;
    res += (1 - 2 * b0 * (xs - x[0])) * a0 * a0 * y[0] + (xs - x[0]) * a0 * a0 * m[0];
    res += (1 - 2 * b1 * (xs - x[1])) * a1 * a1 * y[1] + (xs - x[1]) * a1 * a1 * m[1];

    return res; // 插值基函数与函数值或导数值乘积的线性组合
}
