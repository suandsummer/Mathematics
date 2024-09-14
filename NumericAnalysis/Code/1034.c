/**
 * Problem：用线性元在均匀网格下求解边值问题的数值解
 * Time：2024.05.04
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.141592653898
#define f(x) PI *PI / 2 * sin(PI / 2 * x)
#define N(i, l, r, x) i == 1 ? (r - x) / (r - l) : (x - l) / (r - l)

double Simpson(double left, double right, int flag);                     // Simpson积分公式
void DefKIAndFI(int n, double h, double left, double **KI, double **FI); // 构造单元刚度矩阵KI、单元荷载向量FI
void DefKAndF(int n, double **KI, double **FI, double **K, double *F);   // 构造总刚度矩阵K、总荷载向量F
void swap(int n, double **A, double *b, int r, int k);                   // 交换行列
void Eliminate(int n, double **A, double *b);                            // 列主元消去
void Substitude(int n, double **A, double *b, double *x);                // 回代求解
void Linearmethod(double a, double b, int n, double *xs, double *ys);    // 线性元方法

void OutputMatrix(int m, int n, double **A)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
            printf("%lf ", A[i][j]);
        printf("\n");
    }
}
void OutputVector(int n, double *v)
{
    for (int i = 0; i < n; i++)
    {
        printf("%lf ", v[i]);
        printf("\n");
    }
}

int main()
{
    int N[2];
    double xs[2], ys[2], a = 0, b = 1;

    for (int i = 0; i < 2; ++i)
        scanf("%d,%lf", &N[i], &xs[i]);

    for (int i = 0; i < 2; i++)
        Linearmethod(a, b, N[i], &xs[i], &ys[i]);

    for (int i = 0; i < 2; i++)
        printf("%.4lf\n", ys[i]);

    return 0;
}

double Simpson(double left, double right, int flag)
{
    // Simpson 公式
    double L = right - left;
    double mid = (right + left) / 2;
    double N1 = N(flag, left, right, left);
    double N2 = N(flag, left, right, mid);
    double N3 = N(flag, left, right, right);
    double H = (N3 * f(right) + N1 * f(left) + 4 * N2 * f(mid)) / 6;
    return L * H;
}

void DefKIAndFI(int n, double h, double left, double **KI, double **FI)
{
    double tmp = 1 / h;
    double q = PI * PI / 4;
    KI[0][0] = tmp * (1) + q * h / 3;
    KI[1][1] = tmp * (1) + q * h / 3;
    KI[0][1] = tmp * (-1) + q * h / 6;
    KI[1][0] = tmp * (-1) + q * h / 6;

    for (int i = 0; i < n; i++, left += h)
        for (int j = 0; j < 2; j++)
            FI[i][j] = Simpson(left, left + h, j + 1);
}

void DefKAndF(int n, double **KI, double **FI, double **K, double *F)
{
    for (int i = 0; i <= n; i++)
    {
        F[i] = 0;
        for (int j = 0; j <= n; j++)
            K[i][j] = 0;
    }

    for (int i = 0; i < n; i++)
    {
        K[i][i] += KI[0][0];
        K[i][i + 1] += KI[0][1];
        K[i + 1][i] += KI[1][0];
        K[i + 1][i + 1] += KI[1][1];
        F[i] += FI[i][0];
        F[i + 1] += FI[i][1];
    }
}

void swap(int n, double **A, double *b, int r, int k)
{
    double tmp = b[r];
    b[r] = b[k];
    b[k] = tmp;

    // 交换增广矩阵 (A|b) 之 r,k 两行
    for (int col = 1; col <= n; col++)
    {
        tmp = A[r][col];
        A[r][col] = A[k][col];
        A[k][col] = tmp;
    }
}

void Eliminate(int n, double **A, double *b)
{
    // 消元
    for (int k = 1; k <= n; k++)
    {
        int flag = k;

        // 选择主元
        for (int col = k + 1; col <= n; col++)
            if (fabs(A[col][k]) > fabs(A[flag][k]))
                flag = col;

        // 若主元行与操作行不同，交换两行
        if (k != flag)
            swap(n, A, b, k, flag);

        // 消元
        for (int row = k + 1; row <= n; row++)
        {
            double tmp = A[row][k] / A[k][k];
            A[row][k] = 0;
            for (int coll = k + 1; coll <= n; coll++)
                A[row][coll] -= tmp * A[k][coll];
            b[row] -= b[k] * tmp;
        }
    }
}

void Substitude(int n, double **A, double *b, double *x)
{
    // 回代
    for (int row = n; row > 0; row--)
    {
        b[row] /= A[row][row];
        x[row] = b[row];
        for (int roww = row - 1; roww > 0; roww--)
        {
            b[roww] -= A[roww][row] * x[row];
            A[roww][row] = 0;
        }
    }
}

void Linearmethod(double a, double b, int n, double *xs, double *ys)
{
    // step1：构造单元刚度矩阵 KI、单元荷载向量 FI
    double h = (b - a) / n;
    double **KI = (double **)malloc(sizeof(double *) * 2);
    double **FI = (double **)malloc(sizeof(double *) * n);
    for (int i = 0; i < 2; i++)
        KI[i] = (double *)malloc(sizeof(double) * 2);
    for (int i = 0; i < n; i++)
        FI[i] = (double *)malloc(sizeof(double) * 2);
    DefKIAndFI(n, h, a, KI, FI);

    // printf("\n------KI------\n");
    // OutputMatrix(2, 2, KI);
    // printf("\n------FI------\n");
    // OutputMatrix(n, 2, FI);

    // step2：构造总刚度矩阵 K、总荷载向量 F
    double **K = (double **)malloc(sizeof(double *) * (n + 1));
    double *F = (double *)malloc(sizeof(double) * (n + 1));
    for (int i = 0; i <= n; i++)
        K[i] = (double *)malloc(sizeof(double) * (n + 1));
    DefKAndF(n, KI, FI, K, F);

    // printf("\n------K------\n");
    // OutputMatrix(n + 1, n + 1, K);
    // printf("\n------F------\n");
    // OutputVector(n + 1, F);

    // step3：舍去 0 行 0 列，求解线性方程组 K(1:n)u = F(1:n) （可用追赶法，此处使用列主元消去）
    double *u = (double *)malloc(sizeof(double) * (n + 1));
    Eliminate(n, K, F);
    Substitude(n, K, F, u);
    u[0] = 0;

    // printf("\n------u------\n");
    // OutputVector(n + 1, u);

    // step4：计算 ys
    int flag = 1;
    double left;
    for (left = a; left > *xs || *xs > left + h; left += h)
        flag++;
    *ys = u[flag - 1] * (N(1, left, left + h, *xs)) + u[flag] * (N(2, left, left + h, *xs));

    // printf("\n------flag---left---ys------\n");
    // printf("%d %lf %lf\n", flag, left, *ys);

    // free 回收空间
    for (int i = 0; i < 2; i++)
        free(KI[i]);
    free(KI);
    for (int i = 0; i < n; i++)
        free(FI[i]);
    free(FI);
    for (int i = 0; i <= n; i++)
        free(K[i]);
    free(K);
    free(F);
    free(u);

    return;
}
