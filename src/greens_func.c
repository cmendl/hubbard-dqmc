#include "greens_func.h"
#include "linalg.h"
#include "util.h"
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
	G->mat = (double *)algn_malloc(N*N * sizeof(double));
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
	algn_free(G->mat);
}


//________________________________________________________________________________________________________________________
///
/// \brief Copy Green's function
///
void CopyGreensFunction(const greens_func_t *restrict src, greens_func_t *restrict dst)
{
	assert(src->N == dst->N);

	assume_algned(src->mat);
	assume_algned(dst->mat);
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
	double *Q   = (double *)algn_malloc(N*N * sizeof(double));
	double *tau = (double *)algn_malloc(N   * sizeof(double));   // scalar factors of the elementary reflectors for the matrix Q
	double *d   = (double *)algn_malloc(N   * sizeof(double));
	double *T   = (double *)algn_malloc(N*N * sizeof(double));
	assume_algned(Q);
	assume_algned(tau);
	assume_algned(d);
	assume_algned(T);

	// calculate the imaginary-time flow map
	TimeFlowMap(tsm, slice_shift, Q, tau, d, T);

	// form the matrix D_b^{-1}
	assume_algned(G->mat);
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
	lapack_int *ipiv = algn_malloc(N * sizeof(lapack_int));
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
		G->logdet -= log(t);    // same as log(1/t)

		if (ipiv[i] != i + 1)   // ipiv uses 1-based indices!
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
		G->logdet -= log(t);    // same as log(1/t)
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
	algn_free(ipiv);
	algn_free(T);
	algn_free(d);
	algn_free(tau);
	algn_free(Q);
}


//________________________________________________________________________________________________________________________
///
/// \brief Update the Green's function matrix after a spin flip or phonon update, using the Sherman-Morrison formula
///
void GreenShermanMorrisonUpdate(const double delta, const int N, const int i, double *restrict Gmat)
{
	assume_algned(Gmat);

	int j;

	// copy i-th row of Gmat
	double *c = (double *)algn_malloc(N * sizeof(double));
	assume_algned(c);
	for (j = 0; j < N; j++)
	{
		c[j] = Gmat[i + j*N];
	}
	c[i] -= 1.0;

	// copy i-th column of Gmat
	double *d = (double *)algn_malloc(N * sizeof(double));
	assume_algned(d);
	memcpy(d, &Gmat[N*i], N * sizeof(double));

	// subtract outer Kronecker product from Gmat
	cblas_dger(CblasColMajor, N, N, delta / (1.0 - delta * c[i]), d, 1, c, 1, Gmat, N);

	// clean up
	algn_free(d);
	algn_free(c);
}


//________________________________________________________________________________________________________________________
///
/// \brief Perform the operation B_l G B_l^{-1} with B_l an imaginary-time step generated by the Hamiltonian
///
void GreenTimeSliceWrap(const int N, const double *restrict B, const double *restrict invB, double *restrict Gmat)
{
	// temporary matrix
	double *T = (double *)algn_malloc(N*N * sizeof(double));

	// compute B * G
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, B, N, Gmat, N, 0.0, T, N);

	// compute (B * G) * B^{-1}, store result in Gmat
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, T, N, invB, N, 0.0, Gmat, N);

	// clean up
	algn_free(T);
}


//________________________________________________________________________________________________________________________
///
/// \brief Compute unequal time Green's functions
///
/// \param N            number of lattice sites
/// \param L            number of time steps
/// \param tsm          time step B matrices and precomputed products
/// \param H            temporary array of size L*N x L*N, will be overwritten
/// \param Gtau0        concatenated unequal time Green's functions G(tau,   0) with tau = 0, 1, ..., L-1; array of size L*N x N
/// \param G0tau        concatenated unequal time Green's functions G(0,   tau) with tau = 0, 1, ..., L-1; array of size N x L*N
/// \param Geqlt        concatenated   equal time Green's functions G(tau, tau) with tau = 0, 1, ..., L-1; array of size N x L*N
///
/// Reference:
///   - R. Blankenbecler, D. J. Scalapino, R. L. Sugar\n
///     Monte Carlo calculations of coupled boson-fermion systems. I\n
///     Phys. Rev. D 24, 2278 (1981)
///   - C. Jiang, Z. Bai, R. Scalettar\n
///     A fast selected inversion algorithm for Green's function calculation in many-body quantum Monte Carlo simulations\n
///     IEEE International Parallel and Distributed Processing Symposium, 2016 [http://dx.doi.org/10.1109/IPDPS.2016.69]
///
void ComputeUnequalTimeGreensFunction(const int N, const int L, const time_step_matrices_t *restrict tsm, double *restrict H, double *restrict Gtau0, double *restrict G0tau, double *restrict Geqlt)
{
	int i;

	const int numBprod = tsm->numBprod;

	// set H to the identity matrix
	memset(H, 0, N*N*numBprod*numBprod * sizeof(double));
	for (i = 0; i < numBprod*N; i++)
	{
		H[i + i*numBprod*N] = 1;
	}

	// copy products of B matrices as blocks into H
	for (i = 0; i < numBprod; i++)
	{
		if (i < numBprod-1)
		{
			// lower off-diagonal blocks; note the factor (-1)
			domatcopy(N, N, -1.0, tsm->Bprod[i], N, &H[((i + 1) + numBprod*N*i)*N], numBprod*N);
		}
		else
		{
			// upper right block
			domatcopy(N, N, 1.0, tsm->Bprod[i], N, &H[numBprod*N*(numBprod - 1)*N], numBprod*N);
		}
	}

	// compute the inverse of the block p-cyclic H matrix, with p = numBprod
	BlockCyclicInverse(N, numBprod, H);

	// construct Gtau0 (first column block)
	for (i = 0; i < L; i++)
	{
		if (i % tsm->prodBlen == 0)
		{
			domatcopy(N, N, 1.0, &H[(i/tsm->prodBlen)*N], numBprod*N, &Gtau0[i*N], L*N);
		}
		else
		{
			cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, tsm->B[i-1], N, &Gtau0[(i-1)*N], L*N, 0.0, &Gtau0[i*N], L*N);
		}
	}

	// construct G0tau (first row block)
	for (i = 0; i < L; i++)
	{
		if (i % tsm->prodBlen == 0)
		{
			domatcopy(N, N, 1.0, &H[(i/tsm->prodBlen)*numBprod*N*N], numBprod*N, &G0tau[i*N*N], N);
		}
		else
		{
			if (i == 1)
			{
				// copy first NxN block of H and subtract identity matrix
				double *T = (double *)algn_malloc(N*N * sizeof(double));
				domatcopy(N, N, 1.0, H, numBprod*N, T, N);
				int j;
				for (j = 0; j < N; j++)
				{
					T[j + N*j]--;
				}

				cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, T, N, tsm->invB[i-1], N, 0.0, &G0tau[i*N*N], N);

				algn_free(T);
			}
			else
			{
				cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, &G0tau[(i-1)*N*N], N, tsm->invB[i-1], N, 0.0, &G0tau[i*N*N], N);
			}
		}
	}

	// construct Geqlt (diagonal blocks)
	for (i = 0; i < L; i++)
	{
		if (i % tsm->prodBlen == 0)
		{
			const int ib = (i/tsm->prodBlen);
			domatcopy(N, N, 1.0, &H[(ib + ib*numBprod*N)*N], numBprod*N, &Geqlt[i*N*N], N);
		}
		else
		{
			memcpy(&Geqlt[i*N*N], &Geqlt[(i-1)*N*N], N*N * sizeof(double));
			GreenTimeSliceWrap(N, tsm->B[i-1], tsm->invB[i-1], &Geqlt[i*N*N]);
		}
	}
}
