/**
 * Problem：245 二维卷积操作
 * Time：2024.07.02
 */

#include <stdio.h>
#include <stdlib.h>

// 计算卷积的函数
void convolve2d(int **image, int image_height, int image_width,
                int **kernel, int kernel_height, int kernel_width,
                int **output, int stride, int padding)
{
    int output_height = image_height - kernel_height + 1;
    int output_width = image_width - kernel_width + 1;

    for (int i = 0; i < output_height; ++i)
        for (int j = 0; j < output_width; ++j)
        {
            output[i][j] = 0;
            for (int p = i; p < output_height + i; ++p)
                for (int q = j; q < output_width + j; ++q)
                    output[i][j] += image[p][q] * kernel[p - i][q - j];
        }
}

int main()
{
    // 示例输入图像 (5x5)
    int image_height = 5, image_width = 5;
    int *image_data[5] = {
        (int[]){1, 2, 3, 4, 5},
        (int[]){6, 7, 8, 9, 10},
        (int[]){11, 12, 13, 14, 15},
        (int[]){16, 17, 18, 19, 20},
        (int[]){21, 22, 23, 24, 25}};
    int **image = (int **)image_data;

    // 示例卷积核 (3x3)
    int kernel_height = 3, kernel_width = 3;
    int *kernel_data[3] = {
        (int[]){1, 0, -1},
        (int[]){1, 0, -1},
        (int[]){1, 0, -1}};
    int **kernel = (int **)kernel_data;

    // 计算输出图像的尺寸
    int stride = 1, padding = 0;

    // 初始化输出图像
    int output_height = 3, output_width = 3;
    int **output = (int **)malloc(sizeof(int *) * output_width);
    for (int i = 0; i < output_width; ++i)
        output[i] = (int *)malloc(sizeof(int) * output_height);

    // 执行卷积操作
    convolve2d(image, image_height, image_width, kernel, kernel_height, kernel_width, output, stride, padding);

    // 输出结果
    printf("输入图像：\n");
    for (int i = 0; i < image_height; i++)
    {
        for (int j = 0; j < image_width; j++)
        {
            printf("%2d ", image[i][j]);
        }
        printf("\n");
    }

    printf("\n卷积核：\n");
    for (int i = 0; i < kernel_height; i++)
    {
        for (int j = 0; j < kernel_width; j++)
        {
            printf("%2d ", kernel[i][j]);
        }
        printf("\n");
    }

    printf("\n卷积结果：\n");
    for (int i = 0; i < output_height; i++)
    {
        for (int j = 0; j < output_width; j++)
        {
            printf("%2d ", output[i][j]);
        }
        printf("\n");
    }

    // 释放输出图像的内存
    for (int i = 0; i < output_height; i++)
    {
        free(output[i]);
    }
    free(output);

    return 0;
}
