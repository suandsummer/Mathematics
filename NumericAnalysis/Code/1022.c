/**
 * Problem:1022 梯形公式、Simpson公式和Gauss公式计算积分
 * Time:2024.04.27
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f(double x);                                               // 被积函数
double *compositeGauss(int cells, double left, double right);     // 复合Gauss求积
double Gauss(double left, double right);                          // 对划分的小区间两点求积，两点Gauss公式
double Simpson(double left, double right);                        // 对划分的小区间两点求积，Simpson公式
double trapezoid(double left, double right);                      // 对划分的小区间两点求积，梯形公式
void partition(int cells, double left, double right, double *xs); // 对被积区间做等分
// cells为等分区间的个数

// 其他你需要的函数

int main()
{

    double left, right; // 积分左端点，积分右端点
    double *I;          // 积分值（I[0],I[1],I[2]分别以梯形公式，Simpson公式，两点Gauss公式计算得到）
    int cells = 32;     // 划分的小区间个数

    scanf("%lf%lf", &left, &right);

    I = compositeGauss(cells, left, right);
    I[3] = sin(right) - sin(left);

    // 输出你的答案
    for (int i = 0; i < 4; i++)
        printf("%.10lf\n", I[i]);

    return 0;
}

double f(double x)
{
    return cos(x); // 该实例中， f(x) 取成 cos(x)
}

double *compositeGauss(int cells, double left, double right)
{
    // 返回存放积分值的数组地址，函数声明类型为指针类
    // 利用 malloc 函数为数组 x 动态开辟 cells + 1 个存储空间
    double *x = (double *)malloc(sizeof(double) * (cells + 1));
    static double I[4] = {}; // 存储积分值的数组地址将被返回，应当声明为静态变量，在程序结束时变量消亡

    // 注意：此时若不设定 I[] 为静态变量，在该函数执行结束时该数组变量将被释放，此时将出现返回临时变量地址错误
    // I数组三个分量 I[0],I[1],I[2] 分别存放用梯形公式， Simpson 公式，两点 Gauss 公式计算得到的积分值

    // 对区间做等距划分，将端点信息存入数组 x 中
    // 数组 x 以下标 0 作为起始
    partition(cells, left, right, x);

    // 累加每个小区间的积分值，分别用梯形公式， Simpson 公式，两点 Gauss 公式
    double x1 = left, x2;
    for (int i = 1; i <= cells; i++)
    {
        x2 = x[i];
        I[0] += trapezoid(x1, x2);
        I[1] += Simpson(x1, x2);
        I[2] += Gauss(x1, x2);
        x1 = x2;
    }

    // 利用 free 函数释放数组 x[]
    free(x);

    return I;
}

double Gauss(double left, double right)
{
    // 两点 Gauss 求积公式
    double mid = (left + right) / 2;
    double delta = (right - left) / 2 / sqrt(3);
    double H = (f(mid + delta) + f(mid - delta)) / 2;
    double L = right - left;

    return L * H;
}

double Simpson(double left, double right)
{
    // Simpson 公式
    double mid = (right + left) / 2;
    double H = (f(right) + f(left) + 4 * f(mid)) / 6;
    double L = right - left;
    return L * H;
}

double trapezoid(double left, double right)
{
    // 梯形公式
    double H = (f(right) + f(left)) / 2;
    double L = right - left;
    return L * H;
}

void partition(int cells, double left, double right, double *xs)
{
    // 划分的区间端点以下标 0 作为起始
    // left = x[0] < x[1] < ... < x[cells] = right
    // 计算区间端点
    double h = (right - left) / cells; // 区间长
    for (int i = 0; i <= cells; i++)
        xs[i] = left + h * i;

    return;
}
