/**
 * Problem:1018 Newton插值
 * Time:2024.03.23
 */

#include <stdio.h>
#include <stdlib.h>

void Newton(int n, int m, double x[], double y[], double xs[], double ys[], double d[][25]);

int main()
{

	// 以下程序均舍弃诸数组中的第零个储存空间
	// 若从第零个空间开始记录，请修改 for 循环中初始条件与 malloc 申请内存数量

	int n, m; // n为插值点个数，m为待求点个数
	scanf("%d %d", &n, &m);
	int cells = 5;
	double x[25], y[25], d[25][25]; // x,y 分别记录插值点的横纵坐标；d为 Newton 差商表

	// 利用 malloc 为 xs, ys 数组分配 (m + 1) 个空间
	double *xs = (double *)malloc((m + 1) * sizeof(double));
	double *ys = (double *)malloc((m + 1) * sizeof(double));

	// 依照题意读入各个数据
	for (int i = 1; i <= n; i++)
	{
		scanf("%lf%lf", &x[i], &y[i]);
	}
	for (int i = 1; i <= m; i++)
	{
		scanf("%lf", &xs[i]);
	}

	// Newton 插值
	Newton(n, m, x, y, xs, ys, d);

	for (int i = 1; i <= m; i++)
	{
		printf("%lf %lf\n", xs[i], ys[i]);
	}
	// 释放数组 xs, ys 的空间
	free(xs);
	free(ys);
}

void Newton(int n, int m, double x[], double y[], double xs[], double ys[], double d[][25])
{
	// 对新加入的插值点x[n]更新差商表（行更新）
	for (int i = 1; i <= n; i++)
	{
		d[i][1] = y[i];
		for (int j = 2; j <= i; j++)
		{
			d[i][j] = (d[i][j - 1] - d[i - 1][j - 1]) / (x[i] - x[i - j + 1]);
		}
	}

	// 对每一个待求点 xs[i] 更新计算其 Newton 插值
	for (int i = 1; i <= m; i++)
	{
		ys[i] = d[n][n];
		for (int j = n - 1; j >= 1; j--)
		{
			ys[i] = d[j][j] + (xs[i] - x[j]) * ys[i];
		}
	}

	return;
}
