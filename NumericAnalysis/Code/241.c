/**
 * Problem: 241 FFT and IFFT
 * Time：2024.05.07
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define M_PI 3.14159265358979323846264338327950288

typedef struct
{
    double real;
    double imag;
} complex_t; // 复数

/* Zero has two representations in float point number system
 * It causes problem in OJ system since +0.000000 and -0.000000
 * are different in text.
 * This function convert all values small enought to +0.000000
 * to avoid above problem.
 * It should be mentioned that, in actual practice, the use of
 * this function is typically superfluous.
 * */

double filter_zero(double x);                                                // 判 0 操作
void print_complex_vector(complex_t *data, int n);                           // 输出复数组
complex_t f(double x, double alpha, double beta);                            // 函数 f 产生初始序列
complex_t *new_f_data(double alpha, double beta, double a, double b, int n); // 左端点数组
complex_t complex_mul(complex_t a, complex_t b);                             // 复数乘法
complex_t complex_add(complex_t a, complex_t b);                             // 复数加法
complex_t complex_sub(complex_t a, complex_t b);                             // 复数减法
void fft(complex_t *x, int n);                                               // FFT 快速傅里叶变换
void ifft(complex_t *x, int n);                                              // IFFT 快速傅里叶逆变换

int main()
{
    // NOTE:
    //  Storing 65536-by-65536 complex matrix 'A' requires 64GiB memory.
    //  Computing A*x requires 1.4s by Julia 1.10, 24 threads with 96 cores CPU.
    //  You MUST use FFT to implement DFT.

    // Read parameters:
    int n;
    double alpha, beta;
    scanf("%d", &n); // n=2^k, 1<=k<=16 and k = 16 in some cases.
    scanf("%lf%lf", &alpha, &beta);

    // Generate data:
    complex_t *x = new_f_data(alpha, beta, 0, 2 * M_PI, n);

    // FFT and output result:
    print_complex_vector(x, n);
    fft(x, n);
    print_complex_vector(x, n);
    ifft(x, n);
    for (int i = 0; i < n; ++i)
    {
        x[i].real /= n;
        x[i].imag /= n;
    }
    print_complex_vector(x, n);
    free(x);

    // Generate data:
    x = new_f_data(alpha, beta, 0, 2 * M_PI, n);

    // IFFT and output result:
    print_complex_vector(x, n);
    ifft(x, n);
    for (int i = 0; i < n; ++i)
    {
        x[i].real /= n;
        x[i].imag /= n;
    }
    print_complex_vector(x, n);
    fft(x, n);
    print_complex_vector(x, n);
    free(x);
}

// implement fft and ifft here.

double filter_zero(double x)
{
    if (fabs(x) < 1e-6)
    {
        return fabs(0.0);
    }
    return x;
}

void print_complex_vector(complex_t *data, int n)
{
    int i;
    printf("complex_t[%d] = [\n", n);
    for (i = 0; i < n; i++)
    {
        printf("%.5lf%+.5fim, ", filter_zero(data[i].real), filter_zero(data[i].imag));
        if (i % 4 == 3)
        {
            printf("\n");
        }
    }
    printf("];\n");
}

complex_t f(double x, double alpha, double beta)
{
    // f(x) = cos(ax + b) + sin(ax/2 + b) + i[sin(ax + b) + cos(ax/2 + b)
    complex_t v;
    v.real = cos(alpha * x + beta) + sin(0.5 * alpha * x + beta);
    v.imag = sin(alpha * x + beta) + cos(0.5 * alpha * x + beta);
    return v;
}

complex_t *new_f_data(double alpha, double beta, double a, double b, int n)
{
    // assert(condition); 用于判断条件是否成立，不成立则终止，成立则略过
    assert(a < b);
    assert(n > 0);

    complex_t *data = (complex_t *)malloc(sizeof(complex_t) * n);
    int i;
    double h = (b - a) / n; // not n-1 // 解释：单位根 w0, w1, ... , wn-1
    for (i = 0; i < n; i++)
    {
        data[i] = f(a + i * h, alpha, beta);
    }
    return data;
}

complex_t complex_mul(complex_t a, complex_t b)
{
    complex_t res;
    res.real = a.real * b.real - a.imag * b.imag;
    res.imag = a.real * b.imag + a.imag * b.real;

    return res;
}

complex_t complex_add(complex_t a, complex_t b)
{
    complex_t res;
    res.real = a.real + b.real;
    res.imag = a.imag + b.imag;

    return res;
}

complex_t complex_sub(complex_t a, complex_t b)
{
    complex_t res;
    res.real = a.real - b.real;
    res.imag = a.imag - b.imag;

    return res;
}

void fft(complex_t *x, int n)
{
    // 递归边界
    if (n == 1)
        return;
    int mid = n / 2;
    complex_t x1[mid], x2[mid];

    // 奇偶二分
    for (int i = 0; i < n; i += 2)
    {
        x1[i / 2] = x[i];
        x2[i / 2] = x[i + 1];
    }

    // 分支递归
    fft(x1, mid);
    fft(x2, mid);

    // 合并计算
    complex_t w, t;
    for (int j = 0; j < mid; ++j)
    {
        w.real = cos(2 * M_PI * j / n);
        w.imag = -sin(2 * M_PI * j / n);
        t = complex_mul(w, x2[j]);
        x[j] = complex_add(x1[j], t);
        x[j + mid] = complex_sub(x1[j], t);
    }
}

void ifft(complex_t *x, int n)
{
    // 递归边界
    if (n == 1)
        return;
    int mid = n / 2;
    complex_t x1[mid], x2[mid];

    // 奇偶二分
    for (int i = 0; i < n; i += 2)
    {
        x1[i / 2] = x[i];
        x2[i / 2] = x[i + 1];
    }

    // 分支递归
    ifft(x1, mid);
    ifft(x2, mid);

    // 合并计算
    complex_t w, t;
    for (int j = 0; j < mid; ++j)
    {
        w.real = cos(2 * M_PI * j / n);
        w.imag = sin(2 * M_PI * j / n);
        t = complex_mul(w, x2[j]);
        x[j] = complex_add(x1[j], t);
        x[j + mid] = complex_sub(x1[j], t);
    }
}
