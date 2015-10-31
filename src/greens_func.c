#include "greens_func.h"
#include "linalg.h"
#include <mkl.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>


//________________________________________________________________________________________________________________________
///
/// \brief Allocate memory for the Green's function matrix
///
void AllocateGreensFunction(const int N, greens_func_t *G)
{
	G->mat = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	G->logdet = 0;
	G->sgndet = 1;
	G->N = N;
}


//________________________________________________________________________________________________________________________
///
/// \brief Delete Green's function (free memory)
///
void DeleteGreensFunction(greens_func_t *G)
{
	MKL_free(G->mat);
}


//________________________________________________________________________________________________________________________
///
/// \brief Copy Green's function
///
void CopyGreensFunction(const greens_func_t *restrict src, greens_func_t *restrict dst)
{
	assert(src->N == dst->N);

	__assume_aligned(src->mat, MEM_DATA_ALIGN);
	__assume_aligned(dst->mat, MEM_DATA_ALIGN);
	memcpy(dst->mat, src->mat, src->N*src->N * sizeof(double));

	dst->logdet = src->logdet;
	dst->sgndet = src->sgndet;
}


//________________________________________________________________________________________________________________________
///
/// \brief Construct the Green's function matrix (I + A)^{-1} with A the imaginary-time flow map
///
/// References:
///   - S. R. White, D. J. Scalapino, R. L. Sugar, E. Y. Loh, J. E. Gubernatis, R. T. Scalettar\n
///     Numerical study of the two-dimensional Hubbard model\n
///     Phys. Rev. B 40, 506 (1989)
///   - Z. Bai, C.-R. Lee, R.-C. Li, S. Xu\n
///     Stable solutions of linear systems involving long chain of matrix multiplications\n
///     Linear Algebra Appl. 435, 659-673 (2011)
///
void GreenConstruct(const time_step_matrices_t *restrict tsm, const int slice_shift, greens_func_t *restrict G)
{
	const int N = tsm->N;

	// allocate memory for time flow map
	double *Q   = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	double *tau = (double *)MKL_malloc(N   * sizeof(double), MEM_DATA_ALIGN);	// scalar factors of the elementary reflectors for the matrix Q
	double *d   = (double *)MKL_malloc(N   * sizeof(double), MEM_DATA_ALIGN);
	double *T   = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	__assume_aligned(Q,   MEM_DATA_ALIGN);
	__assume_aligned(tau, MEM_DATA_ALIGN);
	__assume_aligned(d,   MEM_DATA_ALIGN);
	__assume_aligned(T,   MEM_DATA_ALIGN);

	// calculate the imaginary-time flow map
	TimeFlowMap(tsm, slice_shift, Q, tau, d, T);

	// form the matrix D_b^{-1}
	__assume_aligned(G->mat, MEM_DATA_ALIGN);
	memset(G->mat, 0, N*N * sizeof(double));
	int i;
	#pragma ivdep
	for (i = 0; i < N; i++)
	{
		G->mat[i + i*N] = (fabs(d[i]) > 1.0 ? 1.0 / d[i] : 1.0);
	}
	// calculate D_b^{-1} Q^T
	LAPACKE_dormqr(LAPACK_COL_MAJOR, 'R', 'T', N, N, N, Q, N, tau, G->mat, N);

	// calculate D_s T and store result in T
	for (i = 0; i < N; i++)
	{
		if (fabs(d[i]) <= 1.0)
		{
			cblas_dscal(N, d[i], &T[i], N);
		}
	}

	// calculate D_b^{-1} Q^T + D_s T, store result in T
	cblas_daxpy(N*N, 1.0, G->mat, 1, T, 1);

	// perform a LU decomposition of D_b^{-1} Q^T + D_s T
	lapack_int *ipiv = MKL_malloc(N * sizeof(lapack_int), MEM_DATA_ALIGN);
	LAPACKE_dgetrf(LAPACK_COL_MAJOR, N, N, T, N, ipiv);

	// calculate determinant of (D_b^{-1} Q^T + D_s T)^{-1} (D_b^{-1} Q^T)
	G->logdet = 0.0;
	G->sgndet = 1;
	// contribution from (D_b^{-1} Q^T + D_s T)^{-1}
	for (i = 0; i < N; i++)
	{
		double t = T[i + i*N];
		if (t < 0)
		{
			t = -t;
			G->sgndet = -G->sgndet;
		}
		assert(t > 0);
		G->logdet -= log(t);	// same as log(1/t)

		if (ipiv[i] != i + 1)	// ipiv uses 1-based indices!
		{
			G->sgndet = -G->sgndet;
		}

	}
	// contribution from D_b^{-1}
	for (i = 0; i < N; i++)
	{
		if (fabs(d[i]) <= 1.0) {
			continue;
		}

		double t = d[i];
		if (t < 0)
		{
			t = -t;
			G->sgndet = -G->sgndet;
		}
		assert(t > 0);
		G->logdet -= log(t);	// same as log(1/t)
	}
	// Q contributes a factor 1 or -1 to determinant
	for (i = 0; i < N; i++)
	{
		assert(tau[i] >= 0);
		if (tau[i] > 0) {
			G->sgndet = -G->sgndet;
		}
	}

	// calculate (D_b^{-1} Q^T + D_s T)^{-1} (D_b^{-1} Q^T)
	LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', N, N, T, N, ipiv, G->mat, N);

	// clean up
	MKL_free(ipiv);
	MKL_free(T);
	MKL_free(d);
	MKL_free(tau);
	MKL_free(Q);
}


//________________________________________________________________________________________________________________________
///
/// \brief Update the Green's function matrix after a spin flip or phonon update, using the Sherman-Morrison formula
///
void GreenShermanMorrisonUpdate(const double delta, const int N, const int i, double *restrict Gmat)
{
	__assume_aligned(Gmat, MEM_DATA_ALIGN);

	int j;

	// copy i-th row of Gmat
	double *c = (double *)MKL_malloc(N * sizeof(double), MEM_DATA_ALIGN);
	__assume_aligned(c, MEM_DATA_ALIGN);
	for (j = 0; j < N; j++)
	{
		c[j] = Gmat[i + j*N];
	}
	c[i] -= 1.0;

	// copy i-th column of Gmat
	double *d = (double *)MKL_malloc(N * sizeof(double), MEM_DATA_ALIGN);
	__assume_aligned(d, MEM_DATA_ALIGN);
	memcpy(d, &Gmat[N*i], N * sizeof(double));

	// subtract outer Kronecker product from Gmat
	cblas_dger(CblasColMajor, N, N, delta / (1.0 - delta * c[i]), d, 1, c, 1, Gmat, N);

	// clean up
	MKL_free(d);
	MKL_free(c);
}


//________________________________________________________________________________________________________________________
///
/// \brief Perform the operation B_l G B_l^{-1} with B_l an imaginary-time step generated by the Hamiltonian
///
void GreenTimeSliceWrap(const int N, const double *restrict B, const double *restrict invB, double *restrict Gmat)
{
	// temporary matrix
	double *T = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);

	// compute B * G
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, B, N, Gmat, N, 0.0, T, N);

	// compute (B * G) * B^{-1}, store result in Gmat
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, T, N, invB, N, 0.0, Gmat, N);

	// clean up
	MKL_free(T);
}


//________________________________________________________________________________________________________________________
///
/// \brief Compute unequal time Green's functions
///
/// \param N            number of lattice sites
/// \param L            number of time steps
/// \param B            array of B matrices (discrete time step evolution operators), array of length L
/// \param H            temporary array of size L*N x L*N, will be overwritten
/// \param Gtau0        concatenated unequal time Green's functions G(tau,   0) with tau = 0, 1, ..., L-1; array of size L*N x N
/// \param G0tau        concatenated unequal time Green's functions G(0,   tau) with tau = 0, 1, ..., L-1; array of size N x L*N
/// \param Geqlt        concatenated   equal time Green's functions G(tau, tau) with tau = 0, 1, ..., L-1; array of size N x L*N
/// \param Gtau0_avr    G(tau + tau', tau') averaged over tau', with tau, tau' = 0, 1, ..., L-1; array of size L*N x N
/// \param G0tau_avr    G(tau', tau + tau') averaged over tau', with tau, tau' = 0, 1, ..., L-1; array of size N x L*N
/// \param Geqlt_avr    G(tau', tau') averaged over tau', with tau' = 0, 1, ..., L-1; array of size N x N
///
/// Reference:
///   - R. Blankenbecler, D. J. Scalapino, R. L. Sugar\n
///     Monte Carlo calculations of coupled boson-fermion systems. I\n
///     Phys. Rev. D 24, 2278 (1981)
///
void ComputeUnequalTimeGreensFunction(const int N, const int L, const double *const *B, double *restrict H, double *restrict Gtau0, double *restrict G0tau, double *restrict Geqlt, double *restrict Gtau0_avr, double *restrict G0tau_avr, double *restrict Geqlt_avr)
{
	int i, j, k;

	// set H to the identity matrix
	memset(H, 0, N*N*L*L * sizeof(double));
	for (i = 0; i < L*N; i++)
	{
		H[i + i*L*N] = 1;
	}

	// copy B matrices as blocks into H
	for (i = 0; i < L; i++)
	{
		if (i < L-1)
		{
			// lower off-diagonal blocks; note the factor (-1)
			mkl_domatcopy('C', 'N', N, N, -1.0, B[i], N, &H[((i + 1) + L*N*i)*N], L*N);
		}
		else
		{
			// upper right block
			mkl_domatcopy('C', 'N', N, N, 1.0, B[i], N, &H[L*N*(L - 1)*N], L*N);
		}
	}

	// compute the inverse of the block L-cyclic H matrix
	BlockCyclicInverse(N, L, H);

	// Gtau0 is the first column block of H
	memcpy(Gtau0, H, L*N*N * sizeof(double));

	// G0tau is the first row block of H
	mkl_domatcopy('C', 'N', N, L*N, 1.0, H, L*N, G0tau, N);

	// Geqlt contains the diagonal blocks of H
	for (i = 0; i < L; i++)
	{
		mkl_domatcopy('C', 'N', N, N, 1.0, &H[(i + L*N*i)*N], L*N, &Geqlt[i*N*N], N);
	}

	if (Gtau0_avr != NULL)
	{
		// Gtau0_avr is the average of cyclically shifted column blocks of H^{-1}
		memset(Gtau0_avr, 0, L*N*N * sizeof(double));
		for (j = 0; j < L; j++)
		{
			for (i = 0; i < L; i++)
			{
				const int ij = (i + j) % L;

				for (k = 0; k < N; k++)
				{
					cblas_daxpy(N, 1.0, &H[(ij + L*(k + N*j))*N], 1, &Gtau0_avr[(i + k*L)*N], 1);
				}
			}
		}
		// normalize (divide by L)
		cblas_dscal(L*N*N, 1.0/L, Gtau0_avr, 1);
	}

	if (G0tau_avr != NULL)
	{
		// G0tau_avr is the average of cyclically shifted row blocks of H^{-1}
		memset(G0tau_avr, 0, L*N*N * sizeof(double));
		for (j = 0; j < L; j++)
		{
			for (i = 0; i < L; i++)
			{
				const int ij = (i + j) % L;

				for (k = 0; k < N; k++)
				{
					cblas_daxpy(N, 1.0, &H[(ij + L*(k + N*j))*N], 1, &G0tau_avr[(k + j*N)*N], 1);
				}
			}
		}
		// normalize (divide by L)
		cblas_dscal(L*N*N, 1.0/L, G0tau_avr, 1);
	}

	if (Geqlt_avr != NULL)
	{
		// Geqlt_avr is the average of the diagonal blocks of H^{-1}
		memset(Geqlt_avr, 0, N*N * sizeof(double));
		for (i = 0; i < L; i++)
		{
			for (k = 0; k < N; k++)
			{
				cblas_daxpy(N, 1.0, &H[(i + L*(k + N*i))*N], 1, &Geqlt_avr[k*N], 1);
			}
		}
		// normalize (divide by L)
		cblas_dscal(N*N, 1.0/L, Geqlt_avr, 1);
	}
}
