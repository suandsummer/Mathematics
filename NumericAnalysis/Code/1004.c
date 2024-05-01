/**
 * Problem:1004 矩阵范数
 * Time:2024.03.22
 */

#include <stdio.h>
#include <math.h>

double Norm_1(double A[][25], int m, int n);
double Norm_F(double A[][25], int m, int n);
double Norm_inf(double A[][25], int m, int n);

int main()
{

    int n, m;
    double a[25][25];

    scanf("%d%d", &m, &n);

    for (int i = 1; i <= m; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            scanf("%lf", &a[i][j]);
        }
    }

    printf("%lf %lf %lf\n", Norm_1(a, m, n), Norm_F(a, m, n), Norm_inf(a, m, n));
    return 0;
}

double Norm_1(double a[][25], int m, int n)
{
    double max = -1;

    // 计算矩阵 1-范数
    for (int col = 1; col <= n; col++)
    {
        double sum = 0;
        for (int row = 1; row <= m; row++)
        {
            sum += fabs(a[row][col]);
        }
        max = sum > max ? sum : max;
    }

    return max;
}

double Norm_F(double a[][25], int m, int n)
{
    double res = 0;

    // 计算矩阵 Frobenius 范数
    for (int row = 1; row <= m; row++)
    {
        for (int col = 1; col <= n; col++)
        {
            res += pow(a[row][col], 2);
        }
    }
    res = sqrt(res);

    return res;
}

double Norm_inf(double a[][25], int m, int n)
{
    double max = -1;

    // 计算矩阵无穷范数
    for (int row = 1; row <= m; row++)
    {
        double sum = 0;
        for (int col = 1; col <= n; col++)
        {
            sum += fabs(a[row][col]);
        }
        max = sum > max ? sum : max;
    }

    return max;
}
