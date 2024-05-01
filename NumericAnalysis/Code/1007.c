/**
 * Problem:1007 列主元消去法
 * Time:2024.03.22
 */

#include <stdio.h>
#include <math.h>

void swap(int n, double A[][100], double b[], int r, int k);
void Eliminate(int n, double A[][100], double b[]);
void Substitude(int n, double A[][100], double b[], double x[]);

// 以下程序为了和通常所说的矩阵第一行第一列适应，舍去第零行第零列
// 若需要从第零行第零列开始储存，请修改各个循环初始条件与中止条件 (i = 1; ;i <= n) => (i = 0; ; i < n)

void swap1(int n, double A[][100], double b[], int r, int k)
{
    double tmp = b[r];
    b[r] = b[k];
    b[k] = tmp;
    // 交换增广矩阵 (A|b) 之 r,k 两行
    for (int col = 1; col <= n; col++)
    {
        tmp = A[r][col];
        A[r][col] = A[k][col];
        A[k][col] = tmp;
    }
}

void Eliminate1(int n, double A[][100], double b[])
{
    // 消元
    for (int k = 1; k <= n; k++)
    {
        int flag = k;
        // 选择主元
        for (int col = k + 1; col <= n; col++)
            if (fabs(A[col][k]) > fabs(A[flag][k]))
                flag = col;
        // 若主元行与操作行不同，交换两行
        if (k != flag)
            swap1(n, A, b, k, flag);
        // 消元
        for (int row = k + 1; row <= n; row++)
        {
            double tmp = A[row][k] / A[k][k];
            A[row][k] = 0;
            for (int coll = k + 1; coll <= n; coll++)
                A[row][coll] -= tmp * A[k][coll];
            b[row] -= b[k] * tmp;
        }
    }
}

void Substitude1(int n, double A[][100], double b[], double x[])
{
    // 回代
    for (int row = n; row > 0; row--)
    {
        x[row] = b[row] / A[row][row];
        b[row] = 0;
        for (int roww = row - 1; roww > 0; roww--)
        {
            b[roww] -= A[roww][row] * x[row];
            A[roww][row] = 0;
        }
    }
}

int main()
{
    int n;
    double A[100][100] = {}, b[100] = {}, x[100] = {};
    double A1[100][100] = {}, b1[100] = {}, x1[100] = {};

    scanf("%d", &n); // 读入矩阵大小 n , 增广矩阵 (A|b)
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            scanf("%lf", &A[i][j]);
        scanf("%lf", &b[i]);
    }
    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= n; j++)
            A1[i][j] = A[i - 1][j - 1];
        b1[i] = b[i - 1];
    }
    Eliminate1(n, A1, b1);
    Substitude1(n, A1, b1, x1);

    // 输出解向量
    for (int i = 1; i <= n; i++)
        printf("%lf ", x1[i]);

    Eliminate(n, A, b);
    Substitude(n, A, b, x);

    // 输出解向量
    for (int i = 0; i < n; i++)
        printf("%lf ", x[i]);
    return 0;
}

void swap(int n, double A[][100], double b[], int r, int k)
{
    double tmp = b[r];
    b[r] = b[k];
    b[k] = tmp;
    // 交换增广矩阵 (A|b) 之 r,k 两行
    for (int col = 0; col < n; col++)
    {
        tmp = A[r][col];
        A[r][col] = A[k][col];
        A[k][col] = tmp;
    }
}

void Eliminate(int n, double A[][100], double b[])
{
    // 消元
    for (int k = 0; k <= n - 1; k++)
    {
        int flag = k;
        // 选择主元
        for (int col = k + 1; col < n; col++)
            if (fabs(A[col][k]) > fabs(A[flag][k]))
                flag = col;
        // 若主元行与操作行不同，交换两行
        if (k != flag)
            swap(n, A, b, k, flag);
        // 消元
        for (int row = k + 1; row < n; row++)
        {
            double tmp = A[row][k] / A[k][k];
            A[row][k] = 0;
            for (int coll = k + 1; coll < n; coll++)
                A[row][coll] -= tmp * A[k][coll];
            b[row] -= b[k] * tmp;
        }
    }
}

void Substitude(int n, double A[][100], double b[], double x[])
{
    // 回代
    for (int row = n - 1; row > -1; row--)
    {
        x[row] = b[row] / A[row][row];
        b[row] = 0;
        for (int roww = row - 1; roww > -1; roww--)
        {
            b[roww] -= A[roww][row] * x[row];
            A[roww][row] = 0;
        }
    }
}