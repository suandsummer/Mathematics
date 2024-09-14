/**
 * Problem:244 Monte Carlo 法求定积分
 * Time:2024.06.21
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// 定义需要积分的函数，这里以 f(x, y) = x^2 + y^2 为例
double function_to_integrate(double x, double y)
{
    return x * x + y * y;
}

// 蒙特卡洛积分函数（下面这种做法不太标准）
double monte_carlo_integration_2d2(double (*func)(double, double), double x_min, double x_max, double y_min, double y_max, int num_samples)
{
    int cnt = 0;
    double z_max = 2.0;
    double z_min = 0.0;
    for (int i = 0; i < num_samples; ++i)
    {
        double x = (x_max - x_min) * (((double)rand()) / RAND_MAX) + x_min;
        double y = (y_max - y_min) * (((double)rand()) / RAND_MAX) + y_min;
        double z = (z_max - z_min) * (((double)rand()) / RAND_MAX) + z_min;
        if (z < func(x, y))
            cnt++;
    }
    double res = (((double)cnt) / num_samples) * (x_max - x_min) * (y_max - y_min) * (z_max - z_min);
    return res;
}

// 蒙特卡洛积分函数（标准算法）
double monte_carlo_integration_2d(double (*func)(double, double), double x_min, double x_max, double y_min, double y_max, int num_samples)
{
    double sum = 0;
    for (int i = 0; i < num_samples; ++i)
    {
        double x = (x_max - x_min) * (((double)rand()) / RAND_MAX) + x_min;
        double y = (y_max - y_min) * (((double)rand()) / RAND_MAX) + y_min;
        sum += func(x, y);
    }
    return sum / num_samples;
}

int main()
{
    srand(1234);
    int num_samples = 10000;
    double x_min = 0.0, x_max = 1.0;
    double y_min = 0.0, y_max = 1.0;

    // 执行蒙特卡洛积分
    double estimated_integral = monte_carlo_integration_2d(function_to_integrate, x_min, x_max, y_min, y_max, num_samples);

    printf("Estimated integral: %lf\n", estimated_integral);

    return 0;
}
