/**
 * Problem:1010 Jacobi迭代法
 * Time:2024.03.23
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Jocabi(int n, double A[][25], double b[], double x[], double eps, int MaxStep);

int main()
{

    int n;
    double A[25][25], b[25], x[25];
    double eps = 1e-6;
    int MaxStep = 1000;

    scanf("%d", &n);
    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            scanf("%lf", &A[i][j]);
        }
    }
    for (int i = 1; i <= n; i++)
        scanf("%lf", &b[i]);
    for (int i = 1; i <= n; i++)
        scanf("%lf", &x[i]);

    Jocabi(n, A, b, x, eps, MaxStep);
    return 0;
}

void Jocabi(int n, double A[][25], double b[], double x0[], double eps, int MaxStep)
{
    double x[25];
    int cnt = 0;
    double err = eps + 1.0;

    // 中止条件：迭代次数达到 Maxstep 或迭代误差小于 epsilon
    while (err > eps && cnt < MaxStep)
    {
        err = 0;
        // 由 x0[] 计算迭代向量 x[] ，更新迭代向量
        for (int row = 1; row <= n; row++)
        {
            x[row] = b[row];
            for (int col = 1; col <= n; col++)
                if (row != col)
                    x[row] -= A[row][col] * x0[col];
            x[row] /= A[row][row];
            err = err > fabs(x[row] - x0[row]) ? err : fabs(x[row] - x0[row]);
        }
        for (int row = 1; row <= n; row++)
            x0[row] = x[row];
        // 更新迭代次数，记录误差值
        cnt++;
#if 0
        if (cnt == 3)
        {
            printf("\n第%d次迭代值\n", cnt);
            for (int row = 1; row <= n; row++)
                printf("%8.4lf", x0[row]);
            printf("\n");
        }
#endif
    }

    // 判断迭代是否成功，输出结果或错误信息
    if (cnt < MaxStep)
        for (int i = 1; i <= n; i++)
            printf("%lf ", x[i]);
    else
        printf("Fail.");
    return;
}
