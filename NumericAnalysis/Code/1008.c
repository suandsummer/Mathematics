/**
 * Problem:1008 三对角矩阵追赶法
 * Time:2024.03.23
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void read(int a, int b, double arr[]);                                             // 依次读入 arr[a], arr[a + 1], ... , arr[b]
void Doolittle(int n, double d[], double e[], double f[], double l[], double r[]); // 对方阵 A 做 Doolittle 分解
int Chase(int n, double b[], double f[], double l[], double r[], double x[]);      // 追赶过程求解，返回值提示是否求解成功

int main()
{
    int n;
    double b[105], d[105], e[105], f[105], l[105], r[105];
    double x[105];
    int flag;

    scanf("%d", &n);
    read(2, n, d);
    read(1, n, e);
    read(1, n - 1, f);
    read(1, n, b);

    Doolittle(n, d, e, f, l, r);
    flag = Chase(n, b, f, l, r, x);

    if (flag)
    {
        for (int i = 1; i <= n; i++)
        {
            printf("%lf ", x[i]);
        }
    }
    else
    {
        printf("Fail.");
    }

    return 0;
}

void read(int a, int b, double arr[])
{
    for (int i = a; i <= b; i++)
    {
        // 读入 arr[a], arr[a + 1], ... , arr[b]
        scanf("%lf", &arr[i]);
    }
}

void Doolittle(int n, double d[], double e[], double f[], double l[], double r[])
{
    r[1] = e[1];
    for (int i = 2; i <= n; i++)
    {
        // Doolittle分解
        l[i] = d[i] / r[i - 1];
        r[i] = e[i] - l[i] * f[i - 1];
    }
}

int Chase(int n, double b[], double f[], double l[], double r[], double x[])
{
    // 追过程，求解 Ly=b,此处y存在b中
    for (int i = 2; i <= n; i++)
    {
        b[i] -= l[i] * b[i - 1];
    }

    // 赶过程，求解 Ux=y
    if (r[n] == 0)
        return 0;
    x[n] = b[n] / r[n];
    for (int i = n - 1; i >= 1; i--)
    {
        // 若有无法求解的情况，返回 0 值
        x[i] = (b[i] - f[i] * x[i + 1]) / r[i];
    }

    return 1;
}
