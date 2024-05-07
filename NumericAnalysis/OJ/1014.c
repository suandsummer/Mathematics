/**
 * Problem:1014 二分法求根
 * Time:2024.03.22
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double dichotomy(double left, double right, double eps); // 二分法求根函数
double f(double x);

int main()
{

    double left, right; // 给定的左右端点
    double eps = 1e-8;  // 允许误差

    scanf("%lf%lf", &left, &right);

    printf("%lf", dichotomy(left, right, eps));
}

double f(double x)
{
    return sin(x) - tan(x); // 在该实例中， F(x) 取为 sin(x) - tan(x)
}

double dichotomy(double left, double right, double eps)
{

    // 方法一：递归方法求解
    if (fabs(right - left) < eps)
    { // 终止条件：区间长小于 1e-8
        return left;
    }
    double mid = (left + right) / 2;
    if (f(left) * f(mid) < 0)
    {                                     // 利用零点存在定理判断根落在左或右半区间
        return dichotomy(left, mid, eps); // 递归调用二分函数，在左或右半区间继续寻找
    }
    if (f(right) * f(mid) < 0)
    {
        return dichotomy(mid, right, eps);
    }
    /******************************************
        //方法二：循环方式求解
        while(___________){								//满足精度要求时停止循环
            ___________;								//利用零点存在定理判断根在左或右半区间内，更新区间左右端点
        }

        return left;
    ******************************************/
}
