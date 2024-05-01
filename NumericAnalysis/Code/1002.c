/**
 * Problem:1002 Cramer法则
 * Time:2024.03.22
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double determinant(int n, double A[][25]);                 // 求方阵 A[][] 的行列式
int Cramer(int n, double A[][25], double b[], double x[]); // 利用 Cramer 法则求解，返回值代表是否求解成功

// 以下程序为了与通常所说第一行第一列保持一致，舍弃使用第零行第零列
// 若需要从第零行第零列开始储存，请自行修改循环初始与终结条件	(i = 1; ; i <= n) => (i = 0; ; i < n)

int main()
{

    int n; // 矩阵 A 阶数，列向量 b 之维数
    scanf("%d", &n);

    double A[25][25], b[25], x[25];

    for (int i = 1; i <= n; i++)
    { // 读入模块
        for (int j = 1; j <= n; j++)
        {
            scanf("%lf", &A[i][j]);
        }
    }
    for (int i = 1; i <= n; i++)
    {
        scanf("%lf", &b[i]);
    }

    int flag = Cramer(n, A, b, x); // 定义 flag 以记录 Cramer 方法求解成功与否

    if (flag)
    { // 求解成功，输出解向量
        for (int i = 1; i <= n; i++)
        {
            printf("%lf ", x[i]);
        }
    }
    else
    { // 求解失败，抛出错误信息 Fail
        printf("Fail.");
    }
}

double determinant(int n, double A[][25])
{

    if (n == 1)
    {
        return A[1][1]; // 终止条件：一阶方阵，此时行列式值为该元素值
    }

    double det = 0;
    double B[25][25];

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

int Cramer(int n, double A[][25], double *b, double *x)
{

    double D = determinant(n, A); // 矩阵 A 之行列式
    double As[25][25];            // 用向量 b 替换矩阵 A 之各列得到的矩阵

    if (D > 1e-15) // 避免零除错误
    {
        // 先将 A[][] 之各元素赋给 As[][]
        for (int i = 1; i <= n; i++)
        {
            for (int j = 1; j <= n; j++)
            {
                As[i][j] = A[i][j];
            }
        }

        // 计算 x[k]
        for (int k = 1; k <= n; k++)
        {
            // 将矩阵 As 之第 k 列替换为列向量 b
            for (int j = 1; j <= n; j++)
            {
                As[j][k] = b[j];
            }
            // Cramer 法则求解 x[k]
            x[k] = determinant(n, As) / D;
            // 将矩阵 As 还原为
            for (int j = 1; j <= n; j++)
            {
                As[j][k] = A[j][k];
            }
        }

        return 1; // 求解正常，返回真值
    }

    return 0; // 求解失败，返回假值
}
