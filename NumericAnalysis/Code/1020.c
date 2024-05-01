/**
 * Problem:1020 三次样条插值法
 * Time:2024.03.27
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// 给定 n个插值点，计算三次样条插值，并对于 m 个给定的待求点求得其插值函数值
// 使用自然边界条件

void steplen(int n, double *x, double *h);
void cubicspline(int n, int m, double *x, double *y, double *xs, double *ys, double *h);
void ColEliminate(int n, double **A, double *b, double *x);
void swapRow(int n, int r, int k, double **A, double *b);
void compute(int n, int m, double *x, double *y, double *h, double *M, double *xs, double *ys);

int main()
{
    int n, m;
    scanf("%d %d", &n, &m);

    // 使用 malloc 函数为数组开辟空间
    double *x = (double *)malloc(sizeof(double) * n);       // 插值点横坐标
    double *y = (double *)malloc(sizeof(double) * n);       // 插值点纵坐标
    double *xs = (double *)malloc(sizeof(double) * m);      // 待求点横坐标
    double *ys = (double *)malloc(sizeof(double) * m);      // 待求点纵坐标
    double *h = (double *)malloc(sizeof(double) * (n - 1)); // 插值点间距

    for (int i = 0; i < n; i++)
    {
        scanf("%lf %lf", x + i, y + i);
    }
    for (int i = 0; i < m; i++)
    {
        scanf("%lf", xs + i);
    }

    steplen(n, x, h);                   // 计算插值节点间距离
    cubicspline(n, m, x, y, xs, ys, h); // 三次样条插值主函数

    for (int i = 0; i < m; i++)
    {
        printf("%lf %lf\n", xs[i], ys[i]);
    }

    // 使用 free 函数释放各个数组
    free(x);
    free(y);
    free(xs);
    free(ys);
    free(h);

    return 0;
}

void steplen(int n, double *x, double *h)
{
    // 计算插值点间距离
    for (int i = 0; i < n - 1; i++)
        h[i] = x[i + 1] - x[i];
    return;
}

void cubicspline(int n, int m, double *x, double *y, double *xs, double *ys, double *h)
{
    // 使用 malloc 函数为数组开辟空间
    double *a = (double *)malloc(sizeof(double) * n);
    double *b = (double *)malloc(sizeof(double) * n);
    double *c = (double *)malloc(sizeof(double) * n);
    double *d = (double *)malloc(sizeof(double) * n);
    double *M = (double *)malloc(sizeof(double) * n);

    // 计算三对角阵矩阵
    b[0] = 2;
    b[n - 1] = 2;
    c[0] = 1;
    a[n - 1] = 1;
    d[0] = 3 * (y[1] - y[0]) / (x[1] - x[0]); // S''(x_1) = S''(x_n) = 0
    d[n - 1] = 3 * (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);

    for (int i = 1; i < n - 1; i++)
    {
        a[i] = h[i] / (h[i] + h[i - 1]);
        b[i] = 2;
        c[i] = h[i - 1] / (h[i] + h[i - 1]);
        d[i] = 3 * (a[i] * (y[i] - y[i - 1]) / (x[i] - x[i - 1]) + c[i] * (y[i + 1] - y[i]) / (x[i + 1] - x[i]));
    }

    // 列主元消法求解 m（或三对角追赶法求解）

    // 构建矩阵 A
    double **A = (double **)malloc(sizeof(double *) * n);
    for (int i = 0; i < n; i++)
        A[i] = (double *)malloc(sizeof(double) * n);

    for (int i = 0; i < n; i++)
    {
        if (i)
            A[i][i - 1] = a[i];
        A[i][i] = b[i];
        if (i < n - 1)
            A[i][i + 1] = c[i];
    }

    ColEliminate(n, A, d, M);

    // 对于每个待求节点计算插值函数值
    compute(n, m, x, y, h, M, xs, ys);

    for (int i = 0; i < n; i++)
        free(A[i]);
    free(A);
    free(a);
    free(b);
    free(c);
    free(d);
    free(M);

    return;
}

// 列主元消法程序
void ColEliminate(int n, double **A, double *b, double *x)
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
            swapRow(n, k, flag, A, b);
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
void swapRow(int n, int r, int k, double **A, double *b)
{
    double tmp = b[r];
    b[r] = b[k];
    b[k] = tmp;
    for (int i = 0; i < n; i++)
    {
        tmp = A[r][i];
        A[r][i] = A[k][i];
        A[k][i] = tmp;
    }
}

void compute(int n, int m, double *x, double *y, double *h, double *M, double *xs, double *ys)
{
    for (int t = 0; t < m; t++)
    {
        int i;
        for (int j = 0; j < n - 1; j++)
        {
            // 查找 xs[t] 落在何区间内
            if (xs[t] >= x[j] && xs[t] <= x[j + 1])
            {
                i = j;
                break;
            }
        }

        // 计算插值函数值
        double ys1 = pow(xs[t] - x[i + 1], 2) * (h[i] + 2 * (xs[t] - x[i])) * y[i] / pow(h[i], 3);
        double ys2 = pow(xs[t] - x[i], 2) * (h[i] + 2 * (x[i + 1] - xs[t])) * y[i + 1] / pow(h[i], 3);
        double ys3 = pow(xs[t] - x[i + 1], 2) * (xs[t] - x[i]) * M[i] / pow(h[i], 2);
        double ys4 = pow(xs[t] - x[i], 2) * (xs[t] - x[i + 1]) * M[i + 1] / pow(h[i], 2);

        ys[t] = ys1 + ys2 + ys3 + ys4;
    }
}
