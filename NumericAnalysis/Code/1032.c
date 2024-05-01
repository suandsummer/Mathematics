/**
 * Problem:1032 乘幂法
 * Time:2024.04.20
 */

#include <stdio.h>
#include <math.h>
#define eps 1e-6

double EigenValVec(int n, double x0[], double A[][100], int *cnt); // 乘幂法函数
double maxv(int n, double x[]);                                    // 求向量 x 中按模最大者

int main()
{
    int n, cnt = 0;
    double x0[100], A[100][100];

    scanf("%d", &n);
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++)
            scanf("%lf", &A[i][j]);
    for (int i = 1; i <= n; i++)
        scanf("%lf", &x0[i]);

    printf("%lf\n", EigenValVec(n, x0, A, &cnt));

    if (cnt <= 1000)
        for (int i = 1; i <= n; i++)
            printf("%lf ", x0[i]);
    else
        printf("Fail.\n");

    return 0;
}

double EigenValVec(int n, double x0[], double A[][100], int *cnt)
{

    double M = -1, curM = 0; // 记录向量 x 中的按模最大者
    double err = 1;          // 记录迭代误差
    double x[100] = {};      // x = A * x0

    // 迭代步数少于设定的 Maxstep 且迭代误差小于设定的 eps
    while ((*cnt) < 1000 && err > eps)
    {
        // 计算 x
        for (int i = 1; i <= n; ++i)
        {
            x[i] = 0;
            for (int j = 1; j <= n; j++)
                x[i] += A[i][j] * x0[j];
        }

        // 按模最大值归一向量
        curM = maxv(n, x);
        for (int i = 1; i <= n; ++i)
            x[i] /= curM;
        err = fabs(curM - M);
        M = curM;

        // 更新 x0
        for (int i = 1; i <= n; ++i)
        {
            // 更新迭代特征向量与迭代特征值之误差，取最大者
            if (fabs(x0[i] - x[i]) > err)
                err = fabs(x0[i] - x[i]);
            x0[i] = x[i];
        }

        // 记录迭代次数
        (*cnt)++;
    }

    return M; // 返回按模最大特征值
}

double maxv(int n, double x[])
{ // 返回向量 x 中绝对值最大者
    int M = fabs(x[1]), max = 1, cur;
    for (int i = 2; i <= n; i++)
    {
        cur = fabs(x[i]);
        if (M < cur)
        {
            max = i;
            M = cur;
        }
    }
    return x[max];
}
