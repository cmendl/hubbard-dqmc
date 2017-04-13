#include "linalg.h"
#include "dupio.h"
#include <mkl.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>


//________________________________________________________________________________________________________________________
///
/// \brief Compute the matrix exponential of a real symmetric matrix
///
int MatrixExp(const int n, const double *restrict A, double *restrict ret)
{
	__assume_aligned(A, MEM_DATA_ALIGN);

	// eigenvalues
	double *w = MKL_malloc(n * sizeof(double), MEM_DATA_ALIGN);

	// copy 'A' matrix (will be overwritten by eigenvectors)
	double *U = MKL_malloc(n*n * sizeof(double), MEM_DATA_ALIGN);
	__assume_aligned(U, MEM_DATA_ALIGN);
	memcpy(U, A, n*n * sizeof(double));

	int info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', n, U, n, w);
	if (info > 0) {
		duprintf("Intel MKL 'dsyev' failed, error code: %i.\n", info);
		return info;
	}

	// apply exponential to eigenvalues
	int i;
	for (i = 0; i < n; i++)
	{
		w[i] = exp(w[i]);
	}

	// compute U * diag(exp(lambda_i))
	double *UexpD = MKL_malloc(n*n * sizeof(double), MEM_DATA_ALIGN);
	__assume_aligned(UexpD, MEM_DATA_ALIGN);
	memcpy(UexpD, U, n*n * sizeof(double));
	for (i = 0; i < n; i++)
	{
		cblas_dscal(n, w[i], &UexpD[i*n], 1);
	}

	// compute U * diag(exp(lambda_i)) * U^T
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, 1.0, UexpD, n, U, n, 0.0, ret, n);

	// clean up
	MKL_free(UexpD);
	MKL_free(U);
	MKL_free(w);

	return 0;
}


//________________________________________________________________________________________________________________________
///
/// \brief Compute the matrix product A_{num-1} ... A_1 A_0 where each A_i is a n x n matrix
///
void MatrixProductSequence(const int n, const int num, const double *const *restrict A, double *restrict ret)
{
	assert(n > 0);
	assert(num > 0);

	// temporary matrix for calculating products of the A matrices
	double *W = (double *)MKL_malloc(n*n * sizeof(double), MEM_DATA_ALIGN);
	__assume_aligned(W, MEM_DATA_ALIGN);

	// use 'W' and 'ret' as temporary matrices for the products of the A matrices,
	// such that the final results ends up in 'ret'
	double *T1, *T2;
	if ((num & 1) == 0) // if 'num' is even
	{
		T1 = W;
		T2 = ret;
	}
	else
	{
		T1 = ret;
		T2 = W;
	}
	__assume_aligned(T1, MEM_DATA_ALIGN);
	__assume_aligned(T2, MEM_DATA_ALIGN);

	__assume_aligned(A[0], MEM_DATA_ALIGN);
	memcpy(T1, A[0], n*n * sizeof(double));
	int i;
	for (i = 1; i < num; i++)
	{
		// multiply the next A matrix from the left and store in alternating temporary matrices
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, A[i], n, T1, n, 0.0, T2, n);

		// swap T1 and T2 pointers
		{
			double *tmp = T1;
			T1 = T2;
			T2 = tmp;
		}
	}

	// final result should end up in 'ret'
	assert(T1 == ret);

	// clean up
	MKL_free(W);
}


//________________________________________________________________________________________________________________________
///
/// \brief Perform a block p-cyclic QR decomposition
///
/// \param n            block size
/// \param p            cyclic order
/// \param H            block p-cyclic H matrix of size N x N with N = n*p, storing the A_i and B_i blocks; all other entries must be zero at input; H will be overwritten by the block-diagonal entries of Q and R
/// \param tau          scalar factors of the elementary reflectors for the Q blocks at output; vector of length p*n
///
/// Reference:
///     S. Gogolenko, Z. Bai, R. Scalettar\n
///     Structured orthogonal inversion of block p-cyclic matrices on multicore with GPU accelerators\n
///     Euro-Par 2014 Parallel Processing, LNCS 8632, pages 524-535 (2014)\n
///     https://github.com/SGo-Go/BSOFI
///
int BlockCyclicQR(const int n, const int p, double *restrict H, double *restrict tau)
{
	int info;

	const int N  = n*p;

	int i;
	for (i = 0; i < p - 2; i++)
	{
		// perform QR decomposition of [\tilde{A}_i;B_i]; i-th Q block defined by n Householder reflections, but of size 2n x 2n
		info = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, 2*n, n, &H[(1 + N)*i*n], N, tau);
		if (info < 0) {
			duprintf("Intel MKL 'dgeqrf' failed, error code: %i.\n", info);
			return info;
		}

		// multiply by transposed i-th Q block from left
		info = LAPACKE_dormqr(LAPACK_COL_MAJOR, 'L', 'T', 2*n, n, n, &H[(1 + N)*i*n], N, tau, &H[(i + N*(i + 1))*n], N);	if (info < 0) { return info; }
		info = LAPACKE_dormqr(LAPACK_COL_MAJOR, 'L', 'T', 2*n, n, n, &H[(1 + N)*i*n], N, tau, &H[(i + N*(p - 1))*n], N);	if (info < 0) { return info; }

		tau += n;
	}

	// final QR decomposition of the lower right 2n x 2n block
	info = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, 2*n, 2*n, &H[(1 + N)*(p - 2)*n], N, tau);
	if (info < 0) {
		duprintf("Intel MKL 'dgeqrf' failed, error code: %i.\n", info);
		return info;
	}

	return 0;
}


//________________________________________________________________________________________________________________________
///
/// \brief Compute the inverse of an upper triangular matrix with block p-cyclic structure; result is stored in-place
///
/// \param n            block size
/// \param p            cyclic order
/// \param R            block p-cyclic upper triangular matrix of size N x N with N = n*p; will be overwritten by its inverse
///
/// Reference:
///     S. Gogolenko, Z. Bai, R. Scalettar\n
///     Structured orthogonal inversion of block p-cyclic matrices on multicore with GPU accelerators\n
///     Euro-Par 2014 Parallel Processing, LNCS 8632, pages 524-535 (2014)\n
///     https://github.com/SGo-Go/BSOFI
///
int BlockCyclicTriangularInverse(const int n, const int p, double *restrict R)
{
	int info;

	const int N = n*p;

	if (p <= 2)
	{
		info = LAPACKE_dtrtri(LAPACK_COL_MAJOR, 'U', 'N', N, R, N);
		if (info != 0) {
			duprintf("Intel MKL 'dtrtri' failed, error code: %i.\n", info);
			return info;
		}
	}
	else
	{
		// line 2
		info = LAPACKE_dtrtri(LAPACK_COL_MAJOR, 'U', 'N', 3*n, &R[(1 + N)*(p - 3)*n], N);
		if (info != 0) {
			duprintf("Intel MKL 'dtrtri' failed, error code: %i.\n", info);
			return info;
		}

		if (p > 3)
		{
			// multiply last block column by X_{pp}, line 4
			cblas_dtrmm(CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, n*(p - 3), n, 1.0, &R[(1 + N)*(p - 1)*n], N, &R[N*(p - 1)*n], N);

			int i;
			for(i = p - 4; i >= 0; i--)
			{
				// current diagonal block
				double *Rii = &R[(1 + N)*i*n];
				// current off-diagonal block
				double *Rii1 = &R[(i + N*(i + 1))*n];

				// invert diagonal block i, line 3
				info = LAPACKE_dtrtri(LAPACK_COL_MAJOR, 'U', 'N', n, Rii, N);
				if (info != 0) {
					duprintf("Intel MKL 'dtrtri' failed, error code: %i.\n", info);
					return info;
				}

				// line 5(a)
				cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n, n, -1.0, Rii, N, &R[(i + N*(p - 1))*n], N);
				// line 5(b)
				cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n, n, -1.0, Rii, N, Rii1, N);

				// line 7; assuming that blocks not accessed so far have been initialized with zeros
				cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n*(p - i - 2), n, 1.0, Rii1, N, &R[((i + 1) + N*(i + 2))*n], N, 1.0, &R[(i + N*(i + 2))*n], N);
				// line 8
				cblas_dtrmm(CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, n, n, 1.0, &R[(1 + N)*(i + 1)*n], N, Rii1, N);
			}
		}
	}

	return 0;
}


//________________________________________________________________________________________________________________________
///
/// \brief Copy strictly lower triangular matrix entries from A to B, and set these entries in A to zero
///
/// \param m            number of rows in A and B
/// \param n            number of columns in A and B
/// \param A            A matrix
/// \param lda          leading dimension of A, must be >= m
/// \param B            B matrix
/// \param ldb          leading dimension of B, must be >= m
///
static inline void MoveStrictlyLowerTriangularMatrix(const int m, const int n, double *A, const int lda, double *B, const int ldb)
{
	assert(m >= n);

	int j;
	for (j = 0; j < n; j++)
	{
		int i;
		for (i = j + 1; i < m; i++)
		{
			B[i + ldb*j] = A[i + lda*j];
			A[i + lda*j] = 0;
		}
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief Apply transposed orthogonal Q_i blocks from the right to the upper triangular part of the input matrix
///
/// \param n            block size
/// \param p            cyclic order
/// \param A            matrix of size N x N with N = n*p, storing orthogonal Householder reflection data below diagonal
/// \param tau          scalar factors of the elementary reflectors for the Q_i blocks
///
/// Reference:
///     S. Gogolenko, Z. Bai, R. Scalettar\n
///     Structured orthogonal inversion of block p-cyclic matrices on multicore with GPU accelerators\n
///     Euro-Par 2014 Parallel Processing, LNCS 8632, pages 524-535 (2014)\n
///     https://github.com/SGo-Go/BSOFI
///
int BlockCyclicInverseRotation(const int n, const int p, double *restrict A, const double *restrict tau)
{
	const int N = n*p;

	tau += (p - 2)*n;

	// temporary storage for orthogonal Householder reflection data
	// use calloc instead of malloc here because any NaN's on and above the diagonals will be passed
	// into dormqr, and lapack's NaN checking will cause dormqr to immediately return.
	double *Q = (double *)MKL_calloc(4*n*n, sizeof(double), MEM_DATA_ALIGN);

	// copy lower triangular part to 'Q' (orthogonal Householder reflection data) and set entries in 'A' to zero
	MoveStrictlyLowerTriangularMatrix(2*n, 2*n, &A[(1 + N)*(p - 2)*n], N, Q, 2*n);

	// apply 2n x 2n last orthogonal Q_i block from the right
	LAPACKE_dormqr(LAPACK_COL_MAJOR, 'R', 'T', N, 2*n, 2*n, Q, 2*n, tau, &A[N*(p - 2)*n], N);

	int i;
	for (i = p - 3; i >= 0; i--)
	{
		tau -= n;

		// copy lower triangular part of current 2n x n block to 'Q' and set entries in 'A' to zero
		MoveStrictlyLowerTriangularMatrix(2*n, n, &A[(1 + N)*i*n], N, Q, 2*n);

		// apply transposed orthogonal Q_i block from the right
		LAPACKE_dormqr(LAPACK_COL_MAJOR, 'R', 'T', N, 2*n, n, Q, 2*n, tau, &A[N*i*n], N);
	}

	MKL_free(Q);

	return 0;
}


//________________________________________________________________________________________________________________________
///
/// \brief Compute the inverse of the block p-cyclic matrix H
///
/// \param n            block size
/// \param p            cyclic order
/// \param H            block p-cyclic matrix H
///
/// Reference:
///     S. Gogolenko, Z. Bai, R. Scalettar\n
///     Structured orthogonal inversion of block p-cyclic matrices on multicore with GPU accelerators\n
///     Euro-Par 2014 Parallel Processing, LNCS 8632, pages 524-535 (2014)\n
///     https://github.com/SGo-Go/BSOFI
///
int BlockCyclicInverse(const int n, const int p, double *restrict H)
{
	const int N = n*p;

	// allocate memory for scalar factors of the elementary reflectors
	double *tau = MKL_malloc(N * sizeof(double), MEM_DATA_ALIGN);

	int info;
	info = BlockCyclicQR(n, p, H, tau);					if (info != 0) { return info; }
	info = BlockCyclicTriangularInverse(n, p, H);		if (info != 0) { return info; }
	info = BlockCyclicInverseRotation(n, p, H, tau);	if (info != 0) { return info; }

	MKL_free(tau);

	return 0;
}
