/**
 * Problem：实对称矩阵三对角化的 Householder 方法（系统显示超时）
 * Time：2024.06.28
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double norm_inf(int n, double *y);                                          // 计算向量 y 的无穷范数
double norm_2(int n, double *y);                                            // 计算向量 y 的 2 范数
double Householder(int n, double *y, double *v);                            // 计算 Householder 变换, 返回 beta
void Subvector(int n, int k, double **A, double *t);                        // 取向量 t = A(k+1:n,k)
void Vectoru(int n, int k, double **A, double *u, double *v, double gamma); // 计算向量 u
void Vectorw(int m, double *w, double *u, double *v, double gamma);         // 计算向量 w
void Trihouse(int n, double **A);                                           // 三对角化
void Output(int n, double **A);                                             // 输出结果

int main()
{
    double **A;
    int n;

    // malloc 申请内存
    scanf("%d", &n);
    A = (double **)malloc(sizeof(double *) * n);
    for (int i = 0; i < n; ++i)
        A[i] = (double *)malloc(sizeof(double) * n);
    fflush(stdin);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            scanf("%lf", &A[i][j]);

    Trihouse(n, A);

    // 输出结果
    Output(n, A);

    // free 回收内存
    for (int i = 0; i < n; i++)
        free(A[i]);
    free(A);

    return 0;
}
double norm_2(int n, double *y)
{
    double sum = 0;
    for (int i = 0; i < n; i++)
        sum += y[i] * y[i];
    return sqrt(sum);
}

double norm_inf(int n, double *y)
{
    // 计算无穷范数
    double max = fabs(y[0]);
    for (int i = 1; i < n; i++)
        if (max < fabs(y[i]))
            max = fabs(y[i]);

    return max;
}

double Householder(int n, double *y, double *v)
{
    // 计算 Householder 变换, 记录 v[], 返回 beta
    double eta, sigma = 0;

    eta = norm_inf(n, y);
    y[0] /= eta;
    for (int i = 1; i < n; ++i)
    {
        // 避免溢出，做归一化处理
        y[i] /= eta;

        // 计算sigma,v(2:n)
        sigma += y[i] * y[i];
        v[i] = y[i];
    }

    // 计算 beta 和 v[1]
    double beta;
    if (sigma)
    {
        double alpha = sqrt(y[0] * y[0] + sigma);
        if (y[0] <= 0)
            v[0] = y[0] - alpha;
        else
            v[0] = y[0] + alpha;
        beta = 2 * v[0] * v[0] / (sigma + v[0] * v[0]);
        for (int i = 1; i < n; ++i)
            v[i] /= v[0];
        v[0] = 1;
    }
    else
        beta = 0;
    return beta;
}

void Subvector(int n, int k, double **A, double *t)
{
    for (int i = 1; i < n - k; i++)
        t[i - 1] = A[k + i][k];
}

void Vectoru(int n, int k, double **A, double *u, double *v, double gamma)
{
    int m = n - k - 1;
    for (int i = 0; i < m; i++)
    {
        double sum = 0;
        for (int j = 0; j < m; j++)
            sum += A[k + i + 1][k + j + 1] * v[j];
        u[i] = gamma * sum;
    }
}

void Vectorw(int n, double *w, double *u, double *v, double gamma)
{
    double tmp = 0;
    for (int i = 0; i < n; i++)
        tmp += v[i] * u[i];
    tmp *= gamma / 2;
    for (int i = 0; i < n; i++)
        w[i] = u[i] - tmp * v[i];
}

void Trihouse(int n, double **A)
{
    double *u, *w, *t, *v, gamma;
    int flag = 1;

    for (int k = 0; k < n - 2; k++)
    {

        // malloc 申请内存
        int m = n - k - 1;
        u = (double *)malloc(sizeof(double) * m);
        w = (double *)malloc(sizeof(double) * m);
        t = (double *)malloc(sizeof(double) * m);
        v = (double *)malloc(sizeof(double) * m);

        // 三对角化——约化第k行、列
        // printf("\n------第%d次约化------\n", flag++);
        Subvector(n, k, A, t);

        gamma = Householder(m, t, v);

        Vectoru(n, k, A, u, v, gamma);
        Vectorw(n, w, u, v, gamma);

        // printf("\ngamma = %.4lf\nv = ", gamma);
        // for (int i = 0; i < m; i++)
        //     printf("%.4lf ", v[i]);
        // printf("\nu = ");
        // for (int i = 0; i < m; i++)
        //     printf("%.4lf ", u[i]);
        // printf("\nw = ");
        // for (int i = 0; i < m; i++)
        //     printf("%.4lf ", w[i]);

        Subvector(n, k, A, t);
        A[k + 1][k] = -A[k + 1][k] / fabs(A[k + 1][k]) * norm_2(m, t);
        A[k][k + 1] = A[k + 1][k];
        for (int i = k + 2; i < n; i++)
        {
            A[i][k] = 0;
            A[k][i] = 0;
        }
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++)
                A[i + k + 1][j + k + 1] -= v[i] * w[j] + w[i] * v[j];

        // printf("\nA = \n");
        // for (int i = 0; i < n; i++)
        // {
        //     for (int j = 0; j < n; j++)
        //         printf("%.4lf ", A[i][j]);
        //     printf("\n");
        // }

        // free 回收内存
        free(u);
        free(w);
        free(t);
        free(v);
    }
}
void Output(int n, double **A)
{
    // printf("\n------运行结果如下------\n");
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            printf("%lf ", A[i][j]);
        printf("\n");
    }
}

// 4 1 2 1 2 2 2 -1 1 1 -1 1 1 2 1 1 1
// 4 1 2 1 2 2 12.222222 -5.555555 5.555555 1 -5.555555 5.888888 0.111111 2 0.555555 0.111111 1.888888
