/**
 * Problem:1027 最小二乘法
 * Time:2024.04.17
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void leastSquares(int n, int m, double *x, double *y, double *a);                  // 最小二乘法主调用函数
void constructEquation(int n, int m, double *x, double *y, double **A, double *b); // 计算矩阵 A 和列向量 b
void basisFunction(int n, int i, double *x, double *Phi_i);                        // 计算基函数 phi_i 在观测点上的取值做成的列向量 Phi_i
void solveEquation(int m, double **A, double *b, double *a);                       // Doolittle_LU分解方法求解方程组 Aa=b
double dotProduct(int n, double *vec1, double *vec2);                              // 计算向量内积
double calculateAnswer(int m, double *a, double xs);                               // 计算给定点在拟合多项式上的值

int main()
{

    int deg; // 多项式次数
    int n;   // 拟合点个数
    int m;   // 拟合多项式系数个数

    scanf("%d", &n);
    scanf("%d", &deg);

    m = deg + 1;

    double *x; // x[], y[]记录离散点列横纵坐标值
    double *y;
    double *a; // a[] 记录拟合多项式系数

    double xs, ys;

    x = (double *)malloc(sizeof(double) * n); // 利用 malloc 函数为诸数组开辟内存
    y = (double *)malloc(sizeof(double) * n);
    a = (double *)malloc(sizeof(double) * m);

    for (int i = 0; i <= n - 1; i++)
    {
        scanf("%lf %lf", &x[i], &y[i]);
    }
    scanf("%lf", &xs);

    leastSquares(n, m, x, y, a);

    ys = calculateAnswer(m, a, xs);

    for (int i = 0; i <= m - 1; i++)
    {
        printf("%lf ", a[i]);
    }
    printf("\n");
    printf("%lf", ys);

    free(x); // 释放申请的内存
    free(y);
    free(a);
}

void leastSquares(int n, int m, double *x, double *y, double *a)
{

    // 子函数任务一：构造矩阵 A = {(\Phi_i, \Phi_j)} ，构造列向量 b = {(\Phi_i, Y)}
    // 其中， \Phi_i = {\phi_i(x_1), \phi_i(x_2), ..., \phi_i(x_n)}

    // 子函数任务二：求解方程组 Aa = b，则 a 就是所求的最小二乘多项式系数列表

    double **A;
    double *b;

    A = (double **)malloc(sizeof(double *) * m); // 为矩阵 A 和列向量 b 开辟内存
    b = (double *)malloc(sizeof(double) * m);
    for (int i = 0; i < m; i++)
    {
        A[i] = (double *)malloc(sizeof(double) * m);
    }

    constructEquation(n, m, x, y, A, b);

    solveEquation(m, A, b, a);

    for (int i = 0; i < m; i++)
    {
        free(A[i]);
    }
    free(A);
    free(b); // 释放申请的内存
}

void constructEquation(int n, int m, double *x, double *y, double **A, double *b)
{

    // 构造矩阵 A
    double *Phi_i = (double *)malloc(sizeof(double) * n); // 为向量 Phi_i, Phi_j 开辟内存
    double *Phi_j = (double *)malloc(sizeof(double) * n);

    for (int i = 0; i <= m - 1; i++)
    {
        for (int j = 0; j <= m - 1; j++)
        {
            // 计算基函数 \Phi_i, \Phi_j 各个分量
            basisFunction(n, i, x, Phi_i);
            basisFunction(n, j, x, Phi_j);

            // 对基函数向量做点乘，赋值到 A[i][j] 中
            A[i][j] = dotProduct(n, Phi_i, Phi_j);
        }
    }

    // 构造列向量 b
    for (int i = 0; i <= m - 1; i++)
    {
        // 计算基函数向量 \Phi_i 各个分量
        basisFunction(n, i, x, Phi_i);

        // 对 \Phi_i, Y做点乘，赋值到 b[i] 中
        b[i] = dotProduct(n, Phi_i, y);
    }

    free(Phi_i); // 释放申请的内存
    free(Phi_j);
}

void basisFunction(int n, int i, double *x, double *Phi_i)
{

    for (int j = 0; j <= n - 1; j++)
    {

        Phi_i[j] = pow(x[j], i); // 在该实例中，基函数 \phi_i(x) = x^i
    }
}

void solveEquation(int m, double **A, double *b, double *a)
{

    // 利用 Doolittle - LU 分解求解方程组 Aa = b
    double **L;
    double **U;

    // 申请矩阵 L, U 并初始化
    L = (double **)malloc(sizeof(double *) * m);
    U = (double **)malloc(sizeof(double *) * m);
    for (int i = 0; i < m; i++)
    {
        L[i] = (double *)malloc(sizeof(double) * m);
        U[i] = (double *)malloc(sizeof(double) * m);
    }

    // Doolittle分解, A = LU
    for (int i = 0; i < m; i++)
    {
        // 计算矩阵U之第i行第k列元素
        for (int k = i; k < m; k++)
        {
            U[i][k] = A[i][k];
            for (int j = 0; j <= i - 1; j++)
            {
                U[i][k] -= L[i][j] * U[j][k];
            }
        }
        // 计算矩阵L之第i列第k行元素
        for (int k = i; k < m; k++)
        {
            if (k == i)
                L[k][i] = 1;
            else
            {
                L[k][i] = A[k][i];
                for (int j = 0; j <= i - 1; j++)
                {
                    L[k][i] -= L[k][j] * U[j][i];
                }
                L[k][i] /= U[i][i];
            }
        }
    }
#if 0
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            printf("%lf ", A[i][j]);
        }
        printf("\n");
    }
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            printf("%lf ", L[i][j]);
        }
        printf("\n");
    }
    for (int i = 0; i < m; i++)
    {
        for (int j = i; j < m; j++)
        {
            printf("%lf ", U[i][j]);
        }
        printf("\n");
    }
#endif
    // 求解方程组 Aa = LUa = b, 分为 Lz = b, Ua = z两步求解

    double *z = (double *)malloc(sizeof(double) * m);

    for (int col = 0; col < m; col++)
    {
        z[col] = b[col];
        for (int row = col + 1; row < m; row++)
        {
            b[row] -= L[row][col] * z[col];
        }
    }
    for (int col = m - 1; col >= 0; col--)
    {
        a[col] = z[col] / U[col][col];
        for (int row = col - 1; row >= 0; row--)
        {
            z[row] -= U[row][col] * a[col];
        }
    }
    // 释放申请的内存
    for (int i = 0; i < m; i++)
    {
        free(L[i]);
        free(U[i]);
    }
    free(L);
    free(U);
    free(z);
}

double calculateAnswer(int m, double *a, double xs)
{
    double res = a[m - 1]; // 计算 f(xs) = \sum{a_i * x^i}
    for (int i = m - 2; i >= 0; --i)
    {
        res = a[i] + res * xs;
    }
    return res;
}

double dotProduct(int n, double *vec1, double *vec2)
{
    double res = 0;
    for (int i = 0; i < n; i++) // 计算 (vec1,vec2) = \sum{vec1_i * vec_2_i}
    {
        res += vec1[i] * vec2[i];
    }
    return res;
}
