/**
 * Problem:1033 Newton-Cotes求积
 * Time:2024.04.28
 */

#include <stdio.h>
#include <math.h>

double f(double x);                                    // 算例函数
double Newton_Cotes(int n, double left, double right); // Newton-Cotes求积函数
double polynomialintegral(int n, int p[], double x);   // 返回 \int^{x}_{0} p(t) dt
long long factorial(int n);                            // 阶乘函数

int main()
{
    int n;
    double left, right;

    scanf("%d%lf%lf", &n, &left, &right);

    printf("%lf ", Newton_Cotes(n, left, right));

    return 0;
}

double f(double x)
{
    return 1 / (1 + x * x);
}

double Newton_Cotes(int n, double left, double right)
{
    double L = right - left;
    double h = L / n;
    double I = 0;
    for (int i = 0; i <= n; i++)
    {
        // prod 数组储存多项式系数 , prodprime = prod * (t-j)
        int prod[100] = {};
        int prodprime[100] = {};
        prod[0] = 1;

        for (int j = 0; j <= n; j++)
        {
            // 计算多项式 prod
            if (i != j)
            {
                prodprime[0] = prod[0] * (-j);
                for (int k = 1; k <= n; k++)
                    prodprime[k] = prod[k - 1] * 1 + prod[k] * (-j);
                for (int k = 0; k <= n; k++)
                    prod[k] = prodprime[k];
            }
        }

        double int_prod = polynomialintegral(n, prod, n);
        double C = pow(-1, n - i) * int_prod / (n * factorial(i) * factorial(n - i));
        double x = left + h * i;

        I += C * f(x);
    }

    I *= L;

    return I;
}

double polynomialintegral(int n, int p[], double x)
{
    // 对于多项式 p(x) = p[0] + p[1]x + p[2]x^2 + ... + p[n]x^n，计算积分 \int^{x}_{0} p(t) dt
    double sum = 0;
    for (int i = 0; i <= n; ++i)
        sum += p[i] * pow(x, i + 1) / (i + 1);
    return sum;
}

long long factorial(int n)
{
    // 返回阶乘 n!
    long long fac = 1;
    for (int i = 1; i <= n; i++)
        fac *= i;
    return fac;
}
