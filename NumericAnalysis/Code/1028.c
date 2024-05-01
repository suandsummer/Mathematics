/**
 * Problem:1028 Householder变换
 * Time:2024.04.17
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Householder(int n, double y[], double H[][100]); // 计算 Householder 变换矩阵，结果存贮于 H[][] 中
double norm_inf(int n, double y[]);                   // 计算向量 y 的无穷范数

int main()
{

    int n; // 列向量 y 之维数
    double y[100] = {}, H[100][100] = {};

    scanf("%d", &n);
    for (int i = 1; i <= n; i++)
    {
        scanf("%lf", &y[i]);
        H[i][i] = 1;
    }

    Householder(n, y, H);

    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            printf("%lf ", H[i][j]);
        }
        printf("\n");
    }

    return 0;
}

void Householder(int n, double y[], double H[][100])
{
    double eta, sigma = 0, *v;
    v = (double *)malloc(sizeof(double) * (n + 1));

    eta = norm_inf(n, y);
    y[1] /= eta;
    for (int i = 2; i <= n; ++i)
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
        double alpha = sqrt(y[1] * y[1] + sigma);
        if (y[1] <= 0)
            v[1] = y[1] - alpha;
        else
            v[1] = y[1] + alpha;
        beta = 2 * v[1] * v[1] / (sigma + v[1] * v[1]);
        for (int i = 2; i <= n; ++i)
            v[i] /= v[1];
        v[1] = 1;
    }
    else
        beta = 0;

    // 计算变换矩阵 H[][]
    for (int row = 1; row <= n; ++row)
        for (int col = 1; col <= n; ++col)
            H[row][col] -= beta * v[row] * v[col];

    free(v);
}

double norm_inf(int n, double y[])
{
    // 计算列向量 y 的无穷范数
    double max = fabs(y[1]);
    for (int i = 2; i <= n; i++)
        if (max < fabs(y[i]))
            max = fabs(y[i]);
    return max;
}
