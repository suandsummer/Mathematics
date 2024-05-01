/**
 * Problem:1023 Romberg积分法
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f(double x);                                                // 被积函数 f(x)
double trapezoid(double left, double right);                       // 梯形公式
double compositeIntegral(int n, double left, double right);        // 复合 Guass 求积
void partition(int n, double left, double right, double *xs);      // 划分区间端点
void Romberg(int m, int n0, double left, double right, double *I); // Romberg方法求积分

int main()
{
    int m = 4;                                              // 分半加速次数，亦为计算得到积分值个数
    int n0 = 4;                                             // 初始小区间个数
    double left, right;                                     // 积分区间左右端点
    double *I = (double *)malloc(sizeof(double) * (m + 1)); // 利用 malloc 函数声明积分值数组 I[]
    for (int i = 0; i <= m; ++i)
        I[i] = 0;
    scanf("%lf%lf", &left, &right);

    Romberg(m, n0, left, right, I);
    I[m] = sin(right) - sin(left);

    // 输出答案
    for (int i = 0; i <= m; ++i)
        printf("%.10lf\n", I[i]);

    return 0;
}

double f(double x)
{
    return cos(x); // 在该实例中， f(x) 取成 cos(x)
}

void Romberg(int m, int n0, double left, double right, double *I)
{

    double **Q = (double **)malloc(sizeof(double *) * (m + 1)); // 积分表，利用 malloc 函数申请二维数组的列
    for (int i = 0; i <= m; ++i)
        Q[i] = (double *)malloc(sizeof(double) * (m + 1)); // 为列中每一个表头利用 malloc 函数申请行

    // 注意：以下代码从第一行、第一列开始记录，舍弃了第零行和第零列，此时申请的二维数组规模为 (m + 1) * (m + 1)
    // 若从第零行、第零列开始记录数组，需更改下列 for 循环中的循环变量初始条件
    // 若访问了超出申请的地址，将出现指针越界问题而导致运行时错误

    // 对申请的积分表 Q[][] 进行初始化（赋零）
    for (int i = 1; i <= m; ++i)
        for (int j = 1; j <= m; ++j)
            Q[i][j] = 0;
    // 或利用 calloc 函数申请二维数组 Q[][] ，此时数组内元素均为 0，无需额外初始化
    // 利用复合求积公式求得第一次区间划分得到的积分值
    for (int i = 1; i <= m; i++)
        Q[i][1] = compositeIntegral(n0 * pow(2, i - 1) + 1, left, right);

    // 分半加速递推
    for (int j = 2; j <= m; j++)
    { // 循环【列】
        for (int i = j; i <= m; i++)
        { // 循环【行】
            // Romberg递推公式计算 Q[i][j]
            Q[i][j] = (pow(2, 2 * (j - 1)) * Q[i][j - 1] - Q[i - 1][j - 1]) / (pow(2, 2 * (j - 1)) - 1);
        }
    }

    for (int i = 1; i <= m; i++)
    {
        I[i - 1] = Q[i][i]; // 记录积分值到数组 I[]
    }

    // 利用 free 函数释放二维数组 Q[][]
    for (int i = 0; i <= m; ++i)
        free(Q[i]);
    free(Q);

    return;
}

double compositeIntegral(int n, double left, double right)
{ // 参数 n 为划分区间端点的个数

    double *x = (double *)malloc(sizeof(double) * n); // 利用 malloc 函数为数组 x 开辟空间
    double I = 0.0;                                   // 积分值，需初始化为 0

    // 对区间做n等剖分

    partition(n, left, right, x);
    // 利用梯形公式累加小区间积分值
    for (int i = 0; i < n - 1; i++)
        I += trapezoid(x[i], x[i + 1]);

    // 利用 free 函数释放数组 x
    free(x);

    return I;
}

double trapezoid(double left, double right)
{
    // 梯形公式
    double H = (f(right) + f(left)) / 2;
    double L = right - left;
    return L * H;
}

void partition(int n, double left, double right, double *xs)
{
    // 划分得到 n 个端点，记录于 xs 数组中
    // left = x[0] < x[1] < ... < x[n - 1] = right
    // h 为区间步长
    double h = (right - left) / (n - 1);
    // 划分的区间端点以下标 0 作为起始
    for (int i = 0; i <= n - 1; i++)
        xs[i] = left + h * i; // 记录端点信息

    return;
}
