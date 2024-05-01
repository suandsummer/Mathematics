/**
 * Problem:1029 原点位移反幂法
 * Time:2024.04.17
 */

#include <stdio.h>
#include <math.h>

void RevPow(int n, double A[][100], double x[], double k, double eps, double MaxStep); // 反幂法迭代函数
void swap(int n, double A[][100], double b[], int r, int k);                           // 交换行列
void Eliminate(int n, double A[][100], double b[]);                                    // 列主元消去
void Substitude(int n, double A[][100], double b[], double x[]);                       // 回代求解
double max(int n, double y[]);                                                         // 最大值

int main()
{
    int n, MaxStep = 1000;
    double eps = 1e-6, k;
    double A[100][100], x[100];

    scanf("%d", &n);
    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            scanf("%lf", &A[i][j]);
        }
    }
    scanf("%lf", &k);
    for (int i = 1; i <= n; i++)
        scanf("%lf", &x[i]);

    RevPow(n, A, x, k, eps, MaxStep);
    getchar();
    return 0;
}

void RevPow(int n, double A[][100], double x[], double k, double eps, double MaxStep)
{

    int stp = 0;
    double err = 1 + eps, B[100][100], b[100] = {};

    // 原点位移方法，置 A = A - kI
    for (int i = 1; i <= n; i++)
    {
        A[i][i] -= k;
    }

    // 反幂法迭代
    for (int i = 1; i <= n; ++i)
        b[i] = x[i];
    double m = max(n, b);
    while (stp < MaxStep && err > eps)
    {
        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= n; j++)
                B[i][j] = A[i][j];
        Eliminate(n, B, b);
        Substitude(n, B, b, x);
        double cur = max(n, x);
        for (int i = 1; i <= n; ++i)
        {
            x[i] /= cur;
            b[i] = x[i];
        }
        err = fabs(1 / m - 1 / cur);
        m = cur;
        stp++;
    }

    // 输出答案
    if (stp == MaxStep)
        printf("Fail");
    else
    {
        printf("%lf", 1.0 / m + k);
    }
}

void swap(int n, double A[][100], double b[], int r, int k)
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

void Eliminate(int n, double A[][100], double b[])
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
            swap(n, A, b, k, flag);
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

void Substitude(int n, double A[][100], double b[], double x[])
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

double max(int n, double y[])
{
    double max = y[1];
    for (int i = 2; i <= n; i++)
        if (max < y[i])
            max = y[i];
    return max;
}
