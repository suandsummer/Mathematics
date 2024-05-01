/**
 * Problem:1011 Gauss-Seidel方法
 * Time:2024.03.23
 */

#include <stdio.h>
#include <math.h>

void Gauss_Seidel(int n, double a[][100], double b[], double x[]);
// 声明（添加）额外需要的子函数

int main()
{
    int n;
    double a[100][100];
    double b[100], x[100];

    // 读入模块，依次读入 n,A,b,x0
    scanf("%d", &n);
    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            scanf("%lf", &a[i][j]);
        }
    }
    for (int i = 1; i <= n; i++)
        scanf("%lf", &b[i]);
    for (int i = 1; i <= n; i++)
        scanf("%lf", &x[i]);

    Gauss_Seidel(n, a, b, x);
    return 0;
}

void Gauss_Seidel(int n, double a[][100], double b[], double x[])
{
    // 定义你需要的变量
    double eps = 1.0e-6;
    double eps0 = 1 + eps;
    int step = 0, MaxStep = 1000;

    // 判断迭代是否可进行
    double prod = 1;
    for (int i = 1; i <= n; i++)
        if (!(prod *= a[i][i]))
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
                if (row != col)
                    x[row] -= a[row][col] * x[col];
            x[row] = x[row] / a[row][row];
            // 记录误差
            eps0 = eps0 > fabs(x[row] - last) / (1.0 + fabs(x[row])) ? eps0 : fabs(x[row] - last) / (1.0 + fabs(x[row]));
        }
        // 更新迭代次数，更新误差值
        step++;
#if 0
        if (step == 3)
        {
            printf("\n第%d次迭代值\n", step);
            for (int row = 1; row <= n; row++)
                printf("%8.4lf", x[row]);
            printf("\n");
        }
#endif
    }

    // 判断是否迭代成功，输出结果
    // 迭代失败条件：迭代次数到达MaxStep设定值
    if (step < MaxStep)
        for (int row = 1; row <= n; row++)
            printf("%lf ", x[row]);
    else
        printf("Fail.");

    return;
}
