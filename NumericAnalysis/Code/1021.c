/**
 * Problem：B样条插值法
 * Time：2024.05.04
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.1415926535898

double f(double x);
void partition(int n, double left, double right, double *x);                       // 区间等分函数
double Bspline(int n, double left, double right, double *x, double *y, double xs); // B样条插值函数
void swap(int n, double **A, double *b, int r, int k);
void Eliminate(int n, double **A, double *b);
void Substitude(int n, double **A, double *b, double *x);
void ColEliminate(int n, double **A, double *b, double *x); // 列主元消法求解方程组
double Calculate(int n, double *x, double *y, double h, double *M, double left, double xs);
// 计算插值函数值
// 以上 n 为分划得到的小区间个数

int main()
{
    int n = 32;                  // 分划区间个数
    double left = -1, right = 1; // 区间左右端点
    double xs, ys;

    double *x = (double *)malloc(sizeof(double) * (n + 1)); // 利用 malloc 函数为 x,y 数组开辟空间
    double *y = (double *)malloc(sizeof(double) * (n + 1));

    // 读入插值节点
    scanf("%lf", &xs);

    partition(n, left, right, x); // 对区间做 n 等分划

    // 为每个分划区间节点计算函数值
    for (int i = 0; i <= n; i++)
        y[i] = f(x[i]);

    ys = Bspline(n, left, right, x, y, xs);

    // 输出答案，并释放申请的空间
    printf("%lf", ys);
    free(x);
    free(y);

    return 0;
}

double f(double x)
{
    return cos(PI * x); // 在该实例中， f(x) 取为 cos(pi * x)
}

void partition(int n, double left, double right, double *x)
{
    double h = (right - left) / n;

    for (int i = 0; i <= n; ++i)
        x[i] = left + h * i;

    return;
}

double Bspline(int n, double left, double right, double *x, double *y, double xs)
{

    double h = (right - left) / n;                          // 步长
    double *M = (double *)malloc(sizeof(double) * (n + 1)); // 利用 malloc 函数为 M[],d[] 开辟空间
    double *d = (double *)malloc(sizeof(double) * (n + 1));

    double **A = (double **)malloc(sizeof(double *) * (n + 1)); // 利用 malloc 函数为二维数组 A[][] 开辟空间
    for (int i = 0; i <= n; i++)
        A[i] = (double *)malloc(sizeof(double) * (n + 1));

    // 初始化
    double lamda = 0.5, mu = 0.5; // 等距
    for (int i = 0; i <= n; i++)
    {
        M[i] = 0;
        d[i] = 0;
        for (int j = 0; j <= n; j++)
            A[i][j] = 0;
    }

    // 构造系数矩阵 A
    for (int i = 0; i <= n; i++)
    {
        if (i == 0)
        {
            A[0][0] = 1;
            A[0][n] = -1;
        }
        else if (i == 1)
        {
            A[1][1] = 2;
            A[1][2] += lamda;
            A[1][n - 1] += mu;
        }
        else if (i == n)
        {
            A[n][n] = 2;
            A[n][n - 1] += mu;
            A[n][1] += lamda;
        }
        else
        {
            A[i][i] = 2;
            A[i][i - 1] += mu;
            A[i][i + 1] += lamda;
        }
    }

    // 构造列向量 d
    for (int i = 1; i <= n - 1; i++)
        d[i] = 3 * (y[i + 1] - 2 * y[i] + y[i - 1]) / h / h;
    d[n] = 3 * (y[1] - y[0] - y[n] + y[n - 1]) / h / h;

    ColEliminate(n + 1, A, d, M); // 求解 AM = d ，答案存放至数组 M[] 中

    double res = Calculate(n, x, y, h, M, left, xs);

    // 释放申请的空间
    free(M);
    free(d);
    for (int i = 0; i <= n; i++)
        free(A[i]);
    free(A);

    return res; // 返回插值函数值
}

void swap(int n, double **A, double *b, int r, int k)
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

void Eliminate(int n, double **A, double *b)
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

void Substitude(int n, double **A, double *b, double *x)
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

void ColEliminate(int n, double **A, double *b, double *x)
{
    // 列主元消法，自行编写
    Eliminate(n, A, b);
    Substitude(n, A, b, x);
}

double Calculate(int n, double *x, double *y, double h, double *M, double left, double xs)
{
    double ys = 0;
    int i = 0;
    for (int j = 0; j < n; j++)
        if (xs >= x[j] && xs <= x[j + 1])
            i = j;

    // 计算插值函数值 yss
    ys += M[i] * pow(x[i + 1] - xs, 3) / 6 / h;
    ys += M[i + 1] * pow(xs - x[i], 3) / 6 / h;
    ys += (y[i] - M[i] * h * h / 6) * (x[i + 1] - xs) / h;
    ys += (y[i + 1] - M[i + 1] * h * h / 6) * (xs - x[i]) / h;

    return ys;
}
