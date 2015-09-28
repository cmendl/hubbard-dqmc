#include "time_flow.h"
#include "linalg.h"
#include <mkl.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>


//________________________________________________________________________________________________________________________
///
/// \brief Compute a time step B matrix (product of a potential and kinetic energy matrix)
///
inline void ComputeTimeStepMatrix(const kinetic_t *restrict kinetic, const double expV[2], const spin_field_t *restrict s, double *restrict B)
{
	__assume_aligned(B, MEM_DATA_ALIGN);

	const int N = kinetic->N;

	// copy kinetic energy matrix
	__assume_aligned(kinetic->expK, MEM_DATA_ALIGN);
	memcpy(B, kinetic->expK, N*N * sizeof(double));

	// scale the rows by the diagonal potential matrix entries (multiply by diagonal potential matrix from the left)
	int i;
	for (i = 0; i < N; i++)
	{
		cblas_dscal(N, expV[s[i]], &B[i], N);
	}
}

//________________________________________________________________________________________________________________________
///
/// \brief Compute an inverse time step matrix
///
inline void ComputeInverseTimeStepMatrix(const kinetic_t *restrict kinetic, const double expV[2], const spin_field_t *restrict s, double *restrict invB)
{
	__assume_aligned(invB, MEM_DATA_ALIGN);

	const int N = kinetic->N;

	// copy inverse kinetic energy matrix
	__assume_aligned(kinetic->inv_expK, MEM_DATA_ALIGN);
	memcpy(invB, kinetic->inv_expK, N*N * sizeof(double));

	// scale the columns by the inverse diagonal potential matrix entries (multiply by inverse diagonal potential matrix from the right)
	int i;
	for (i = 0; i < N; i++)
	{
		cblas_dscal(N, expV[1-s[i]], &invB[i*N], 1);
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief Compute a time step B matrix (product of a potential and kinetic energy matrix), including phonons 
///
inline void ComputePhononTimeStepMatrix(const kinetic_t *restrict kinetic, const double expV[2], const spin_field_t *restrict s, const double *restrict expX, double *restrict B)
{
	__assume_aligned(B, MEM_DATA_ALIGN);

	const int N = kinetic->N;

	// copy kinetic energy matrix
	__assume_aligned(kinetic->expK, MEM_DATA_ALIGN);
	memcpy(B, kinetic->expK, N*N * sizeof(double));

	// scale the rows by the diagonal potential matrix entries (multiply by diagonal potential matrix from the left)
	int i;
	for (i = 0; i < N; i++)
	{
		cblas_dscal(N, expV[s[i]] * expX[i], &B[i], N);
	}
}

//________________________________________________________________________________________________________________________
///
/// \brief Compute an inverse time step matrix including phonons
///
inline void ComputeInversePhononTimeStepMatrix(const kinetic_t *restrict kinetic, const double expV[2], const spin_field_t *restrict s, const double *restrict expX, double *restrict invB)
{
	__assume_aligned(invB, MEM_DATA_ALIGN);

	const int N = kinetic->N;

	// copy inverse kinetic energy matrix
	__assume_aligned(kinetic->inv_expK, MEM_DATA_ALIGN);
	memcpy(invB, kinetic->inv_expK, N*N * sizeof(double));

	// scale the columns by the inverse diagonal potential matrix entries (multiply by inverse diagonal potential matrix from the right)
	int i;
	for (i = 0; i < N; i++)
	{
		cblas_dscal(N, expV[1-s[i]] / expX[i], &invB[i*N], 1);
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief Allocate memory for the time step B matrices and precomputed products
///
void AllocateTimeStepMatrices(const int N, const int L, const int prodBlen, time_step_matrices_t *restrict tsm)
{
	assert(L % prodBlen == 0);

	tsm->N = N;
	tsm->L = L;
	tsm->prodBlen = prodBlen;
	tsm->numBprod = L / prodBlen;
	assert(tsm->numBprod > 0);

	tsm->B    = (double **)MKL_malloc(L * sizeof(double *), MEM_DATA_ALIGN);
	tsm->invB = (double **)MKL_malloc(L * sizeof(double *), MEM_DATA_ALIGN);

	tsm->Bprod = (double **)MKL_malloc(tsm->numBprod * sizeof(double *), MEM_DATA_ALIGN);

	// ensure that each individual matrix is aligned at a 'MEM_DATA_ALIGN' boundary
	int l;
	for (l = 0; l < L; l++)
	{
		tsm->B[l]    = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
		tsm->invB[l] = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	}
	for (l = 0; l < tsm->numBprod; l++)
	{
		tsm->Bprod[l] = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief Free memory of the time step matrices structure
///
void DeleteTimeStepMatrices(time_step_matrices_t *restrict tsm)
{
	int l;
	for (l = 0; l < tsm->numBprod; l++)
	{
		MKL_free(tsm->Bprod[l]);
	}
	for (l = 0; l < tsm->L; l++)
	{
		MKL_free(tsm->invB[l]);
		MKL_free(tsm->B[l]);
	}

	MKL_free(tsm->Bprod);
	MKL_free(tsm->invB);
	MKL_free(tsm->B);

	tsm->L = 0;
}


//________________________________________________________________________________________________________________________
///
/// \brief Copy time step matrices, assuming that memory for the target structure has already been allocated
///
void CopyTimeStepMatrices(time_step_matrices_t *restrict src, time_step_matrices_t *restrict dst)
{
	// all dimensions must agree
	assert(dst->N ==  src->N);
	assert(dst->L ==  src->L);
	assert(dst->prodBlen == src->prodBlen);
	assert(dst->numBprod == src->numBprod);

	const int N = src->N;

	const double *srcM;
	double *dstM;

	int l;
	for (l = 0; l < src->L; l++)
	{
		srcM = src->B[l];
		dstM = dst->B[l];
		__assume_aligned(srcM, MEM_DATA_ALIGN);
		__assume_aligned(dstM, MEM_DATA_ALIGN);
		memcpy(dstM, srcM, N*N * sizeof(double));

		srcM = src->invB[l];
		dstM = dst->invB[l];
		__assume_aligned(srcM, MEM_DATA_ALIGN);
		__assume_aligned(dstM, MEM_DATA_ALIGN);
		memcpy(dstM, srcM, N*N * sizeof(double));
	}
	for (l = 0; l < src->numBprod; l++)
	{
		srcM = src->Bprod[l];
		dstM = dst->Bprod[l];
		__assume_aligned(srcM, MEM_DATA_ALIGN);
		__assume_aligned(dstM, MEM_DATA_ALIGN);
		memcpy(dstM, srcM, N*N * sizeof(double));
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief Initialize time step B matrices and products of subsequent 'tsm->prodBlen' matrices
///
void InitTimeStepMatrices(const kinetic_t *restrict kinetic, const double expV[2], const spin_field_t *restrict s, time_step_matrices_t *restrict tsm)
{
	int l;
	const int N = kinetic->N;
	assert(tsm->N == N);

	// calculate B and B^{-1} matrices
	for (l = 0; l < tsm->L; l++)
	{
		ComputeTimeStepMatrix       (kinetic, expV, &s[l*N], tsm->B[l]);
		ComputeInverseTimeStepMatrix(kinetic, expV, &s[l*N], tsm->invB[l]);
	}

	// form products of the B matrices
	for (l = 0; l < tsm->numBprod; l++)
	{
		MatrixProductSequence(N, tsm->prodBlen, (const double **)(tsm->B + l*tsm->prodBlen), tsm->Bprod[l]);
	}
}

//________________________________________________________________________________________________________________________
///
/// \brief Initialize time step B matrices and products of subsequent 'tsm->prodBlen' matrices, taking phonons into account
///
void InitPhononTimeStepMatrices(const kinetic_t *restrict kinetic, const double expV[2], const spin_field_t *restrict s, const double *restrict expX, time_step_matrices_t *restrict tsm)
{
	int l;
	const int N = kinetic->N;
	assert(tsm->N == N);

	// calculate B and B^{-1} matrices
	for (l = 0; l < tsm->L; l++)
	{
		ComputePhononTimeStepMatrix       (kinetic, expV, &s[l*N], &expX[l*N], tsm->B[l]);
		ComputeInversePhononTimeStepMatrix(kinetic, expV, &s[l*N], &expX[l*N], tsm->invB[l]);
	}

	// form products of the B matrices
	for (l = 0; l < tsm->numBprod; l++)
	{
		MatrixProductSequence(N, tsm->prodBlen, (const double **)(tsm->B + l*tsm->prodBlen), tsm->Bprod[l]);
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief Update time step B matrices for time slice 'l', and re-compute products of subsequent matrices
/// if 'l' refers to the last B matrix of a product interval
///
void UpdateTimeStepMatrices(const kinetic_t *restrict kinetic, const double expV[2], const spin_field_t *restrict s, const int l, time_step_matrices_t *restrict tsm)
{
	assert(0 <= l && l < tsm->L);

	ComputeTimeStepMatrix       (kinetic, expV, &s[l*tsm->N], tsm->B[l]);
	ComputeInverseTimeStepMatrix(kinetic, expV, &s[l*tsm->N], tsm->invB[l]);

	// re-compute products of B matrices if 'l' refers to the last B matrix of a product interval
	if (l % tsm->prodBlen == tsm->prodBlen - 1)
	{
		const int lstart = l - (tsm->prodBlen - 1);
		MatrixProductSequence(tsm->N, tsm->prodBlen, (const double **)(tsm->B + lstart), tsm->Bprod[lstart / tsm->prodBlen]);
	}
}

//________________________________________________________________________________________________________________________
///
/// \brief Update time step B matrices for time slice 'l' taking phonons into account,
/// and re-compute products of subsequent matrices if 'l' refers to the last B matrix of a product interval
///
void UpdatePhononTimeStepMatrices(const kinetic_t *restrict kinetic, const double expV[2], const spin_field_t *restrict s, const double *restrict expX, const int l, time_step_matrices_t *restrict tsm)
{
	assert(0 <= l && l < tsm->L);

	ComputePhononTimeStepMatrix       (kinetic, expV, &s[l*tsm->N], &expX[l*tsm->N], tsm->B[l]);
	ComputeInversePhononTimeStepMatrix(kinetic, expV, &s[l*tsm->N], &expX[l*tsm->N], tsm->invB[l]);

	// re-compute products of B matrices if 'l' refers to the last B matrix of a product interval
	if (l % tsm->prodBlen == tsm->prodBlen - 1)
	{
		const int lstart = l - (tsm->prodBlen - 1);
		MatrixProductSequence(tsm->N, tsm->prodBlen, (const double **)(tsm->B + lstart), tsm->Bprod[lstart / tsm->prodBlen]);
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief Calculate the imaginary-time flow map generated by the Hamiltonian using Lie-Trotter splitting;
/// the output is represented as Q * diag(d) * T with Q orthogonal
///
/// \param tsm                  time step matrices structure
/// \param slice_shift          cyclic time slice shift of the B matrices; must be a multiple of 'tsm->prodBlen'
/// \param Q                    orthogonal matrix stored in the format used by the LAPACK 'dgeqp3' function
/// \param tau                  scalar factors of the elementary reflectors for the matrix Q; refer to the LAPACK 'dgeqp3' function
/// \param d                    entries of the diagonal matrix in the output representation
/// \param T                    T matrix in the output representation
///
/// References:
///   - S. R. White, D. J. Scalapino, R. L. Sugar, E. Y. Loh, J. E. Gubernatis, R. T. Scalettar\n
///     Numerical study of the two-dimensional Hubbard model\n
///     Phys. Rev. B 40, 506 (1989)
///   - Z. Bai, C.-R. Lee, R.-C. Li, S. Xu\n
///     Stable solutions of linear systems involving long chain of matrix multiplications\n
///     Linear Algebra Appl. 435, 659-673 (2011)
///
void TimeFlowMap(const time_step_matrices_t *restrict tsm, const int slice_shift, double *restrict Q, double *restrict tau, double *restrict d, double *restrict T)
{
	__assume_aligned(Q,   MEM_DATA_ALIGN);
	__assume_aligned(tau, MEM_DATA_ALIGN);
	__assume_aligned(d,   MEM_DATA_ALIGN);
	__assume_aligned(T,   MEM_DATA_ALIGN);

	int i;
	const int N = tsm->N;

	assert(tsm->L == tsm->prodBlen * tsm->numBprod);

	// divide 'slice_shift' by product length of the B matrices
	assert(slice_shift >= 0);
	assert(slice_shift % tsm->prodBlen == 0);	// must be a multiple of 'tsm->prodBlen'
	const int prod_slice_shift = slice_shift / tsm->prodBlen;

	// temporary matrix for the QR decompositions of the product of subsequent 'tsm->prodBlen' B matrices
	double *W = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);

	// use input 'Q' as temporary matrix for storing the orthogonal rotations
	double *U = Q;

	if ((tsm->numBprod & 1) == 1)	// if 'numBprod' is odd
	{
		// swap W and U pointers to ensure that the final orthogonal rotation ends up in 'Q'
		double *tmp = W;
		W = U;
		U = tmp;
	}

	__assume_aligned(W, MEM_DATA_ALIGN);
	__assume_aligned(U, MEM_DATA_ALIGN);

	// permutation array for QR decomposition with pivoting
	lapack_int* jpvt = (lapack_int *)MKL_malloc(N * sizeof(lapack_int), MEM_DATA_ALIGN);
	__assume_aligned(jpvt, MEM_DATA_ALIGN);

	// temporary array storing the inverse entries of 'd' or columns of 'T'
	double *v = (double *)MKL_malloc(N * sizeof(double), MEM_DATA_ALIGN);
	__assume_aligned(v, MEM_DATA_ALIGN);

	// initial QR decomposition
	{
		// copy precomputed product of subsequent B matrices
		const double *srcBprod = tsm->Bprod[prod_slice_shift % tsm->numBprod];
		__assume_aligned(srcBprod, MEM_DATA_ALIGN);
		memcpy(W, srcBprod, N*N * sizeof(double));

		// perform a QR-decomposition with column pivoting of the product of subsequent B matrices
		memset(jpvt, 0, N * sizeof(lapack_int));	// all columns of the input matrix are "free" columns
		LAPACKE_dgeqp3(LAPACK_COL_MAJOR, N, N, W, N, jpvt, tau);

		// extract diagonal entries
		for (i = 0; i < N; i++)
		{
			d[i] = W[i + i*N];
			// avoid division by zero
			if (d[i] == 0) {
				d[i] = 1;
			}
			v[i] = 1.0 / d[i];
		}

		// construct the matrix T_1
		memset(T, 0, N*N * sizeof(double));
		int j;
		for (j = 0; j < N; j++)
		{
			// indices start at zero in C
			jpvt[j]--;

			for (i = 0; i <= j; i++)
			{
				T[i + jpvt[j]*N] = v[i] * W[i + j*N];
			}
		}
	}

	int l;
	for (l = 1; l < tsm->numBprod; l++)
	{
		// swap W and U pointers
		{
			double *tmp = W;
			W = U;
			U = tmp;
		}
		__assume_aligned(W, MEM_DATA_ALIGN);
		__assume_aligned(U, MEM_DATA_ALIGN);

		// copy product of the next 'tsm->prodBlen' B matrices
		const double *srcBprod = tsm->Bprod[(l + prod_slice_shift) % tsm->numBprod];
		__assume_aligned(srcBprod, MEM_DATA_ALIGN);
		memcpy(W, srcBprod, N*N * sizeof(double));

		// multiply by previous orthogonal matrix from the right, overwriting the 'W' matrix
		LAPACKE_dormqr(LAPACK_COL_MAJOR, 'R', 'N', N, N, N, U, N, tau, W, N);
		// scale by previous diagonal matrix
		for (i = 0; i < N; i++)
		{
			cblas_dscal(N, d[i], &W[i*N], 1);
		}

		// perform a QR-decomposition of 'W' with column pivoting
		memset(jpvt, 0, N * sizeof(lapack_int));	// all columns of the input matrix are "free" columns
		LAPACKE_dgeqp3(LAPACK_COL_MAJOR, N, N, W, N, jpvt, tau);

		// adjust pivoting indices
		for (i = 0; i < N; i++)
		{
			// indices start at zero in C
			jpvt[i]--;
		}

		// extract diagonal entries
		#pragma ivdep
		for (i = 0; i < N; i++)
		{
			d[i] = W[i + i*N];
			// avoid division by zero
			if (d[i] == 0) {
				d[i] = 1;
			}
			v[i] = 1.0 / d[i];
		}

		// multiply inverse diagonal matrix from the left with the upper triangular R matrix
		int j;
		for (j = 0; j < N; j++)
		{
			for (i = 0; i <= j; i++)
			{
				W[i + j*N] *= v[i];
			}
		}

		// multiply current T_l matrix with the product of the previous T_l's:
		// first, apply pivoting permutation to the rows of the product of the previous T_l's
		for (j = 0; j < N; j++)
		{
			// store permuted column in 'v'
			for (i = 0; i < N; i++)
			{
				v[i] = T[jpvt[i] + j*N];
			}
			// copy 'v' back into 'T'
			memcpy(&T[j*N], v, N*sizeof(double));
		}
		// second, multiply with upper triangular matrix, overwriting 'T' in-place;
		// we do not assume that upper triangular matrix is unit triangular since entries might be zero
		cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, N, N, 1.0, W, N, T, N);
	}

	// pointers should match; 'W' still contains the final orthogonal rotation information in its lower triangular part
	assert(W == Q);

	// clean up
	MKL_free(v);
	MKL_free(jpvt);
	MKL_free(U);	// 'U' points to the memory location of the original 'W' matrix
}
