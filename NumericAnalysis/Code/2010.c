/**
 * Problem:2010 Jacobi迭代
 * Time:2024.03.23
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{
    double A[3][3] = {8, -1, 1, 2, 10, -1, 1, 1, -5};
    double b[3] = {1, 4, 3};
    double x[3], x0[3] = {0.125, 0.4, -0.6};
    double eps = 1e-6;
    double err = 1 + eps;
    int cnt = 0, Maxstep = 1000;

    while (err > eps && cnt < Maxstep)
    {
        err = 0;
        for (int row = 0; row < 3; row++)
        {
            x[row] = b[row];
            for (int col = 0; col < 3; col++)
                if (row != col)
                    x[row] -= A[row][col] * x0[col];
            x[row] /= A[row][row];
            err = err > fabs(x[row] - x0[row]) ? err : fabs(x[row] - x0[row]);
        }
        for (int row = 0; row < 3; row++)
            x0[row] = x[row];
        cnt++;
    }

    if (cnt < Maxstep)
    {
        printf("%d\n", cnt);
        for (int row; row < 3; row++)
            printf("%lf ", x[row]);
    }
    else
        printf("Fali.\n");

    return 0;
}