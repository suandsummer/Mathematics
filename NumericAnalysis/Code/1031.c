/**
 * Problem:1031 正交变换法
 * Time:2024.04.17
 */

#include <stdio.h>
#include <math.h>

void QR(int m, int n, double A[][100], double b[], double c[]);                 // 对矩阵 A 做 QR 分解, 计算列向量 c = Q'b
double Householder(int n, double y[], double v[]);                              // 计算 Householder 变换, 返回 beta
double norm_inf(int n, double y[]);                                             // 计算向量 y 的无穷范数
void Solve(int m, double A[][100], double c[], double x[]);                     // 求解方程组 Rx = c
void MatrixMultiply(int m, int n, double A[][100], double B[][100]);            // 矩阵乘法，结果存于B
void HouseholderMatrix(int j, int m, double beta, double v[], double H[][100]); // 计算H矩阵

int main()
{
    int m, n; // 矩阵 A[][] ,m 行 n 列
    double A[100][100], b[100], c[100] = {};
    double x[100];

    scanf("%d %d", &m, &n);
    for (int i = 1; i <= m; i++)
        for (int j = 1; j <= n; j++)
            scanf("%lf", &A[i][j]);
    for (int i = 1; i <= m; i++)
        scanf("%lf", &b[i]);

    QR(m, n, A, b, c);
    Solve(n, A, c, x);

    for (int i = 1; i <= n; i++)
        printf("%lf ", x[i]);
    return 0;
}

void QR(int m, int n, double A[][100], double b[], double c[])
{
    double Q[100][100] = {}, v[100] = {}, beta;

    // Q 由累乘各步 Householder 矩阵得到，初值置为 I
    for (int i = 1; i <= m; i++)
        Q[i][i] = 1;

    for (int j = 1; j <= n; j++)
    {
        if (j < m)
        {
            // Householder变换计算 v[], beta
            double a[100] = {};

            // a = A(j:m, j)
            for (int i = j; i <= m; i++)
                a[i - j + 1] = A[i][j];
            beta = Householder(m - j + 1, a, v);

            // H = I - beta * v * v', 即第 j 步得到的 Householder 矩阵
            double H[100][100] = {};
            // 计算变换矩阵 H
            HouseholderMatrix(j, m, beta, v, H);

            // 计算 A(j:m, j:n)
            MatrixMultiply(m, n, H, A);

            // 计算第 j 步得到的正交矩阵 Q[][], 不再计算 d[]
            MatrixMultiply(m, m, H, Q);
        }
    }

    // 计算 c = Q_1' * b, Q_1 = Q(1:m, 1:m)
    for (int i = 1; i <= m; i++)
    {
        c[i] = 0;
        for (int k = 1; k <= m; k++)
        {
            c[i] += Q[i][k] * b[k];
        }
    }
}

void Solve(int n, double A[][100], double c[], double x[])
{
    // 求解上三角方程组 Rx = c
    for (int j = n; j >= 1; j--)
    {
        x[j] = c[j] / A[j][j];
        for (int i = j - 1; i >= 1; i--)
        {
            c[i] -= A[i][j] * x[j];
        }
    }
}

double Householder(int n, double y[], double v[])
{

    // 计算 Householder 变换, 记录 v[], 返回 beta
    double eta, sigma = 0;

    eta = norm_inf(n, y);
    y[1] /= eta;
    for (int i = 2; i <= n; ++i)
    {
        // 避免溢出，做归一化处理
        y[i] /= eta;

        // 计算sigma,v(2:n)
        sigma += y[i] * y[i];
        v[i] = y[i];
    }

    // 计算 beta 和 v[1]
    double beta;
    if (sigma)
    {
        double alpha = sqrt(y[1] * y[1] + sigma);
        if (y[1] <= 0)
            v[1] = y[1] - alpha;
        else
            v[1] = y[1] + alpha;
        beta = 2 * v[1] * v[1] / (sigma + v[1] * v[1]);
        for (int i = 2; i <= n; ++i)
            v[i] /= v[1];
        v[1] = 1;
    }
    else
        beta = 0;
    return beta;
}

double norm_inf(int n, double y[])
{
    // 计算无穷范数
    double max = fabs(y[1]);
    for (int i = 2; i <= n; i++)
        if (max < fabs(y[i]))
            max = fabs(y[i]);
    return max;
}

void MatrixMultiply(int m, int n, double A[][100], double B[][100])
{
    double T[m + 1][n + 1];
    for (int i = 1; i <= m; i++)
    {
        for (int k = 1; k <= n; k++)
        {
            T[i][k] = 0;
            for (int l = 1; l <= m; l++)
            {
                T[i][k] += A[i][l] * B[l][k];
            }
        }
    }
    for (int i = 1; i <= m; i++)
        for (int j = 1; j <= n; j++)
            B[i][j] = T[i][j];
}

void HouseholderMatrix(int j, int m, double beta, double v[], double H[][100])
{
    for (int i = 1; i <= m; i++)
        H[i][i] = 1;
    for (int row = j; row <= m; ++row)
        for (int col = j; col <= m; ++col)
            H[row][col] -= beta * v[row - j + 1] * v[col - j + 1];
}