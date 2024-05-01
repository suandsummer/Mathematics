/**
 * Problem:1017 Lagrange插值
 * Time:2024.03.23
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Lagrange(int n, int m, double x[], double y[], double xs[], double ys[]);

int main()
{

    int n, m;
    scanf("%d %d", &n, &m);

    // 利用malloc函数动态开辟堆内存
    double *x = (double *)malloc(sizeof(double) * (n + 1));
    double *y = (double *)malloc(sizeof(double) * (n + 1));
    double *xs = (double *)malloc(sizeof(double) * (m + 1));
    double *ys = (double *)malloc(sizeof(double) * (m + 1));

    // 读入插值节点 x_i,y_i，读入待求点 xs_i
    for (int i = 1; i <= n; i++)
    {
        scanf("%lf%lf", &x[i], &y[i]);
    }
    for (int i = 1; i <= m; i++)
    {
        scanf("%lf", &xs[i]);
    }

    Lagrange(n, m, x, y, xs, ys);

    // 输出插值结果
    for (int i = 1; i <= m; i++)
        printf("%lf %lf\n", xs[i], ys[i]);

    // 释放申请内存
    free(x);
    free(y);
    free(xs);
    free(ys);

    return 0;
}

void Lagrange(int n, int m, double x[], double y[], double xs[], double ys[])
{

    for (int k = 1; k <= m; k++)
    {
        // k = 1...m ，求 ys[k] 的值
        for (int i = 1; i <= n; i++)
        {
            // 计算插值基函数 l_i(x)
            double l = 1.0;
            for (int j = 1; j <= n; j++)
            {
                if (j != i)
                {
                    l *= (xs[k] - x[j]) / (x[i] - x[j]);
                }
            }
            // 计算第 i 个插值节点对 ys[k] 的贡献
            ys[k] += y[i] * l;
        }
    }
}
