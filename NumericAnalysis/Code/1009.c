/**
 * Problem:1009 一般形式迭代法
 * Time:2024.03.23
 */

#include <stdio.h>
#include <math.h>

void iteration(int n, double H[][100], double g[], double x[], double x0[], double eps, int MaxStep);

int main()
{

    int n, MaxStep = 1000;
    double eps = 1e-6;
    double x0[100], x[100]; // x = H * x0 + g
    double H[100][100], g[100];

    scanf("%d", &n);
    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            scanf("%lf", &H[i][j]);
        }
    }
    for (int i = 1; i <= n; i++)
    {
        scanf("%lf", &g[i]);
    }
    for (int i = 1; i <= n; i++)
    {
        scanf("%lf", &x0[i]);
    }

    iteration(n, H, g, x, x0, eps, MaxStep);

    return 0;
}

void iteration(int n, double H[][100], double g[], double x[], double x0[], double eps, int MaxStep)
{

    // 添加你需要定义的变量
    double eps0 = 1 + eps;
    int step = 0;
    for (int row = 1; row <= n; row++)
        x[row] = x0[row];
    // 迭代
    while (step < MaxStep && eps0 > eps)
    { // 迭代进行的条件：迭代次数不超过MaxStep,误差不小于1e-6
        // 计算x = H * x0 + g
        eps0 = 0;
        for (int row = 1; row <= n; row++)
        {
            double sum = 0;
            for (int col = 1; col <= n; col++)
            {
                sum += H[row][col] * x0[col];
            }
            x[row] = sum + g[row];
            // 更新误差值，记录迭代次数，更新x0
            eps0 = (eps0 > x[row] - x0[row]) ? eps0 : x[row] - x0[row];
        }
        for (int row = 1; row <= n; row++)
            x0[row] = x[row];
        step++;
    }

    // 判断是否迭代成功，输出结果
    // 迭代失败条件：迭代次数到达MaxStep设定值
    if (step < MaxStep)
    {
        for (int row = 1; row <= n; row++)
        {
            printf("%lf ", x[row]);
        }
    }
    else
    {
        printf("Fail.");
    }

    return;
}
