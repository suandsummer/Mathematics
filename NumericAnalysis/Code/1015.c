/**
 * Problem:1015 Newton迭代法
 * Time:2024.03.22
 */

#include <stdio.h>
#include <math.h>
#define E 0.000001

double iterate(double x, double a, double b, double c, double d); // 迭代函数

#define E 0.000001

double iterate(double x, double a, double b, double c, double d); // 迭代函数

double f(double x, double a, double b, double c, double d)
{
    return a * pow(x, 3) + b * pow(x, 2) + c * x + d;
}

double ff(double x, double a, double b, double c)
{
    return 3 * a * pow(x, 2) + 2 * b * x + c;
}
// 此处声明（添加）需要的子函数

int main()
{

    double x, a, b, c, d;

    scanf("%lf %lf %lf %lf %lf", &x, &a, &b, &c, &d);
    printf("%lf", iterate(x, a, b, c, d));

    return 0;
}

double iterate(double x, double a, double b, double c, double d)
{
    double x2 = x;
    double x1, f1, f2;
    do
    {
        f1 = f(x2, a, b, c, d);
        x1 = x2;
        x2 = x2 - f1 / ff(x2, a, b, c);
        f2 = f(x2, a, b, c, d);
    } while (fabs(f1 - f2) > E);
    return x2;
}