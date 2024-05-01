/**
 * Problem:1006 Crout分解
 * Time:2024.03.22
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Crout(int n, double A[][100], double L[][100], double U[][100]);

// 以下代码舍弃了第零行第零列，若从第零行第零列开始记录需要更改 for 循环内初值条件

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

    Crout(n, A, L, U);

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

double determinant(int n, double A[][100])
{

    if (n == 1)
    {
        return A[1][1]; // 终止条件：一阶方阵，此时行列式值为该元素值
    }

    double det = 0;
    double B[100][100] = {};

    for (int i = 1; i <= n; i++)
    { // 按第一行展开，对每个 A[1][i] 计算其代数余子式
        for (int j = 1; j <= n; j++)
        {              // 循环 A[][] 之各列，记录余子式 B[][]
            if (j < i) // 循环到的第 j 列在第 i 列左侧
            {
                for (int k = 1; k <= n - 1; k++)
                {
                    B[k][j] = A[k + 1][j];
                }
            }
            if (j > i) // 循环到的第 j 列在第 i 列右侧
            {
                for (int k = 1; k <= n - 1; k++)
                {
                    B[k][j - 1] = A[k + 1][j];
                }
            }
        }
        det += pow(-1, i + 1) * A[1][i] * determinant(n - 1, B); // 累加代数余子式与该元素之积
    }

    return det;
}

void Judge(int n, double A[][100])
{
    if (!n)
        return;
    if (determinant(n, A))
        Judge(n - 1, A);
    else
    {
        printf("Fail");
        exit(1);
    }
}

void Crout(int n, double A[][100], double L[][100], double U[][100])
{
    // 判断 Crout 分解是否能够进行
    Judge(n, A);

    // 进行Crout分解
    for (int i = 1; i <= n; i++)
    {
        // 计算矩阵 L
        for (int k = i; k <= n; k++)
        {
            L[k][i] = A[k][i];
            for (int j = 1; j <= i - 1; j++)
            {
                L[k][i] -= L[k][j] * U[j][i];
            }
        }
        // 计算矩阵 U
        for (int k = i; k <= n; k++)
        {
            if (k == i)
                U[k][i] = 1;
            else
            {
                U[i][k] = A[i][k];
                for (int j = 1; j <= i - 1; j++)
                {
                    U[i][k] -= L[i][j] * U[j][k];
                }
                U[i][k] /= L[i][i];
            }
        }
    }
    return;
}
