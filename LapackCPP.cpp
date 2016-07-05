#include <mkl_types.h>
#include <mkl_lapack.h>
#include"mkl.h"
//定义解双精度复共轭矩阵的本征值问题的Lapack库的C++接口
void heev_cpp(MKL_Complex16 *Hamilton, double EigenValue[], int N)
{
	int info;
	int lda;
	int lwork;
	char jobz;
	char uplo;
	double rwork[3 * N - 2];
	MKL_Complex16 work[2 * N - 1];

	lda = N;
	lwork = 2 * N - 1;
	jobz = 'V';
	uplo = 'L';

	ZHEEV(&jobz, &uplo, &N, Hamilton, &lda, EigenValue, work, &lwork, rwork, &info);
}

void inverse(MKL_Complex16* A, int N)
{
	int m = N;
	int n = N;
	int lda = N;
	int ipiv[N];
	int info;
	info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, m, n, A, lda, ipiv);
	info = LAPACKE_zgetri(LAPACK_ROW_MAJOR, m, A, lda, ipiv);

}
