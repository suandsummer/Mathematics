/**
 * Problem:1005 Doolittle分解
 * Time:2024.3.22
 */

#include <stdio.h>
#include <math.h>

void Doolittle(int n, double A[][100], double L[][100], double U[][100]);

int main()
{
    int n;
    double A[100][100] = {}, L[100][100] = {}, U[100][100] = {};

    scanf("%d", &n);
    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            scanf("%lf", &A[i][j]);
        }
    }
    Doolittle(n, A, L, U);

    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            printf("%lf ", L[i][j]);
        }
        printf("\n");
    }
    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            printf("%lf ", U[i][j]);
        }
        printf("\n");
    }
    getchar();
    return 0;
}
void Doolittle(int n, double A[][100], double L[][100], double U[][100])
{
    // Doolittle分解
    for (int i = 1; i <= n; i++)
    {
        // 计算矩阵U之第i行第k列元素
        for (int k = i; k <= n; k++)
        {
            U[i][k] = A[i][k];
            for (int j = 1; j <= i - 1; j++)
            {
                U[i][k] -= L[i][j] * U[j][k];
            }
        }
        // 计算矩阵L之第i列第k行元素
        for (int k = i; k <= n; k++)
        {
            if (k == i)
                L[k][i] = 1;
            else
            {
                L[k][i] = A[k][i];
                for (int j = 1; j <= i - 1; j++)
                {
                    L[k][i] -= L[k][j] * U[j][i];
                }
                L[k][i] /= U[i][i];
            }
        }
    }
    return;
}