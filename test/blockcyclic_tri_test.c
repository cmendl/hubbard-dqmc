#include "linalg.h"
#include "random.h"
#include "util.h"
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>


int BlockCyclicTriTest()
{
	const int n = 4;
	const int p = 5;

	const int N = p*n;

	// artificial UNIX time
	time_t itime = 1420000241;

	// random generator seed; multiplicative constant from Pierre L'Ecuyer's paper
	randseed_t seed;
	Random_SeedInit(1865811235122147685LL * (uint64_t)itime, &seed);

	// allocate R matrix
	double *R = (double *)algn_malloc(N*N * sizeof(double));
	memset(R, 0, N*N * sizeof(double));

	int i, j, k;

	// define random non-zero blocks of R matrix
	for (i = 0; i < p; i++)
	{
		// reference to current R_ii block
		double *Rii = &R[(i + N*i)*n];

		// diagonal blocks are upper triangular
		for (k = 0; k < n; k++)
		{
			for (j = 0; j <= k; j++)
			{
				Rii[j + N*k] = Random_GetUniform(&seed) - 0.5;
			}
		}
	}
	for (i = 0; i < p - 1; i++)
	{
		// reference to current R_{i,i+1} block
		double *Rii1 = &R[(i + N*(i + 1))*n];

		for (k = 0; k < n; k++)
		{
			for (j = 0; j < n; j++)
			{
				Rii1[j + N*k] = Random_GetUniform(&seed) - 0.5;
			}
		}
	}
	for (i = 0; i < p - 2; i++)
	{
		// reference to current R_{i,p-1} block
		double *Rip = &R[(i + N*(p - 1))*n];

		for (k = 0; k < n; k++)
		{
			for (j = 0; j < n; j++)
			{
				Rip[j + N*k] = Random_GetUniform(&seed) - 0.5;
			}
		}
	}

	// copy R for later reference
	double *R_ref = (double *)algn_malloc(N*N * sizeof(double));
	memcpy(R_ref, R, N*N * sizeof(double));

	// compute the inverse of R
	printf("Computing the inverse of an upper triangular matrix with block p-cyclic structure...\n");
	BlockCyclicTriangularInverse(n, p, R);
	//LAPACKE_dtrtri(LAPACK_COL_MAJOR, 'U', 'N', N, R, N);

	// calculate matrix product Rinv*R
	double *M = (double *)algn_malloc(N*N * sizeof(double));
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, R, N, R_ref, N, 0.0, M, N);

	// subtract identity matrix
	for (i = 0; i < N; i++)
	{
		M[(1 + N)*i]--;
	}

	// entrywise absolute error (M should be zero)
	double err = 0;
	for (i = 0; i < N*N; i++)
	{
		err = fmax(err, fabs(M[i]));
	}
	printf("Largest entrywise absolute error: %g\n", err);

	// clean up
	algn_free(M);
	algn_free(R_ref);
	algn_free(R);

	return (err < 8e-14 ? 0 : 1);
}
