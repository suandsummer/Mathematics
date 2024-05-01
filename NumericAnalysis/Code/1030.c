/**
 * Problem:1030 正则化方法
 * Time:2024.04.17
 */

#include <stdio.h>
#include <math.h>

void product(int n, int m, double A[][100], double C[][100], double d[], double b[]);
void Cholesky(int n, double A[][100], double L[][100]);
void Chase(int m, double L[][100], double d[], double x[], double y[]);

int main()
{

    int n, m;
    double A[100][100], C[100][100], d[100], b[100], L[100][100];
    double x[100], y[100];

    scanf("%d %d", &n, &m); // 矩阵 A[][] n 行 m 列
    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= m; j++)
        {
            scanf("%lf", &A[i][j]);
        }
    }
    for (int i = 1; i <= n; i++)
    {
        scanf("%lf", &b[i]);
    }

    product(n, m, A, C, d, b);
    Cholesky(m, C, L);
    Chase(m, L, d, x, y);

    for (int i = 1; i <= m; i++)
    {
        printf("%lf ", x[i]);
    }

    return 0;
}

void product(int n, int m, double A[][100], double C[][100], double d[], double b[])
{
    // 计算 C = A' * A
    for (int i = 1; i <= m; i++)
    {
        for (int j = 1; j <= m; j++)
        {
            C[i][j] = 0;
            for (int k = 1; k <= n; k++)
                C[i][j] += A[k][i] * A[k][j];
        }
    }
    // 计算 d = A' * b
    for (int i = 1; i <= m; i++)
    {
        d[i] = 0;
        for (int j = 1; j <= n; j++)
        {
            d[i] += A[j][i] * b[j];
        }
    }
}

void Cholesky(int n, double A[][100], double L[][100])
{
    // Cholesky 分解正定对称阵 A
    for (int j = 1; j <= n; j++)
    {
        for (int i = j; i <= n; i++)
        {
            if (i == j)
            {
                L[i][i] = A[i][i];
                for (int k = 1; k <= i - 1; k++)
                {
                    L[i][i] -= L[i][k] * L[i][k];
                }
                L[i][i] = sqrt(L[i][i]);
            }
            else
            {
                L[i][j] = A[i][j];
                for (int k = 1; k <= j - 1; k++)
                {
                    L[i][j] = L[i][k] * L[j][k];
                }
                L[i][j] /= L[j][j];
            }
        }
    }
}

void Chase(int n, double L[][100], double d[], double x[], double y[])
{
    // 求解方程组 Ly = d
    for (int j = 1; j <= n; j++)
    {
        y[j] = d[j] / L[j][j];
        for (int i = j + 1; i <= n; i++)
        {
            d[i] -= L[i][j] * y[j];
        }
    }

    // 求解方程组 L'x = y
    for (int j = n; j >= 1; j--)
    {
        x[j] = y[j] / L[j][j];
        for (int i = j - 1; i >= 1; i--)
        {
            y[i] -= L[j][i] * x[j];
        }
    }
}
