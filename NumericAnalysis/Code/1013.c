/**
 * Problem:1013 SOR方法
 * Time:2024.03.27
 */

#include <stdio.h>
#include <math.h>

// 额外需要的子函数
void SOR(int n, double A[][100], double b[], double x[], double w);

int main()
{

    int n;
    double A[100][100], b[100], x[100];
    double w;

    // 读入模块
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
    scanf("%lf", &w);

    SOR(n, A, b, x, w);

    return 0;
}

void SOR(int n, double A[][100], double b[], double x[], double w)
{
    // 需要定义的变量
    double eps = 1.0e-6;
    double eps0 = 1 + eps;
    int step = 0, MaxStep = 1000;

    // 判断SOR方法是否可进行
    if (w <= 0 || w >= 2)
    {
        printf("Fail.");
        return;
    }
    double prod = 1;
    for (int i = 1; i <= n; i++)
        if (!(prod *= A[i][i]))
        {
            printf("Fail.");
            return;
        }

    // 迭代进行的条件：迭代次数不超过MaxStep,误差不小于1e-6
    while (eps0 > eps && step < MaxStep)
    {
        // 构造滚动迭代序列(Seidel技巧)
        eps0 = 0;
        for (int row = 1; row <= n; row++)
        {
            double last = x[row];
            x[row] = b[row];
            for (int col = 1; col <= n; col++)
            {
                if (row == col)
                    x[row] -= A[row][col] * last;
                else
                    x[row] -= A[row][col] * x[col];
            }
            x[row] = last + w * x[row] / A[row][row];
            // 记录误差
            if (1)
                eps0 = eps0 > fabs(x[row] - last) / (1.0 + fabs(x[row])) ? eps0 : fabs(x[row] - last) / (1.0 + fabs(x[row]));
            else
                eps0 = eps0 > fabs(x[row] - last) ? eps0 : fabs(x[row] - last);
        }
        // 更新迭代次数，更新误差值
        step++;
    }

    // 判断是否迭代成功，输出结果
    // 迭代失败条件：迭代次数到达MaxStep设定值
    if (step < MaxStep)
    {
        for (int row = 1; row <= n; row++)
            printf("%lf ", x[row]);
        printf("\n%d", step);
    }
    else
        printf("Fail.");

    return;
}