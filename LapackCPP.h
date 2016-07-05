#ifndef _CPP_LAPACK95_H_
#define _CPP_LAPACK95_H_

#include <mkl_types.h>

//计算双精度复数型矩阵的本征值问题
void heev_cpp(MKL_Complex16 *Hamilton, double EigenValue[], int N);

//求逆矩阵，类型待定
void inverse(MKL_Complex16* A, int N);

#endif
