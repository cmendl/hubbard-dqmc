#include "linalg.h"
#include "random.h"
#include <mkl.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>


int BlockCyclicInvTest()
{
	const int n = 4;
	const int p = 5;

	const int N = p*n;

	// artificial UNIX time
	time_t itime = 1420000241;

	// random generator seed; multiplicative constant from Pierre L'Ecuyer's paper
	randseed_t seed;
	Random_SeedInit(1865811235122147685LL * (uint64_t)itime, &seed);

	// allocate H matrix
	double *H = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	memset(H, 0, N*N * sizeof(double));

	int i, j, k;

	// define random A and B blocks
	for (i = 0; i < p; i++)
	{
		// reference to current A block
		double *A = &H[(i + N*i)*n];

		// reference to current B block
		double *B;
		if (i < p - 1) {
			B = &H[((i + 1) + N*i)*n];
		}
		else {
			B = &H[N*(p - 1)*n];
		}

		for (k = 0; k < n; k++)
		{
			for (j = 0; j < n; j++)
			{
				A[j + N*k] = Random_GetUniform(&seed) - 0.5;
				B[j + N*k] = Random_GetUniform(&seed) - 0.5;
			}
		}
	}

	// copy H for later reference
	double *H_ref = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	memcpy(H_ref, H, N*N * sizeof(double));

	// Computing the inverse of the block p-cyclic H matrix
	printf("Compute the inverse of a block p-cyclic matrix...\n");
	BlockCyclicInverse(n, p, H);

	// calculate matrix product Hinv*H
	double *M = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, H, N, H_ref, N, 0.0, M, N);

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
	MKL_free(M);
	MKL_free(H_ref);
	MKL_free(H);

	return (err < 4e-15 ? 0 : 1);
}
