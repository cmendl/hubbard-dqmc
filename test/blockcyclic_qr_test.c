#include "linalg.h"
#include "random.h"
#include "util.h"
#include <mkl.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>


int BlockCyclicQRTest()
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

	// allocate tau vector
	double *tau = (double *)MKL_malloc(N * sizeof(double), MEM_DATA_ALIGN);;

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

	// perform block p-cyclic QR decomposition
	printf("Performing a block p-cyclic QR decomposition...\n");
	BlockCyclicQR(n, p, H, tau);

	// copy entries of 'H' to 'Q' for preserving orthogonal Q blocks
	double *Q = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	memcpy(Q, H, N*N * sizeof(double));

	// set sub-diagonal entries of R (overwritten entries of H) to zero
	for (i = 0; i < N - 1; i++)
	{
		memset(&H[i + 1 + N*i], 0, (N - 1 - i)*sizeof(double));
	}

	// calculate matrix product Q*R and store result back in H
	LAPACKE_dormqr(LAPACK_COL_MAJOR, 'L', 'N', 2*n, N, 2*n, &Q[((p - 2) + N*(p - 2))*n], N, tau + (p - 2)*n, &H[(p - 2)*n], N);
	for (i = p - 3; i >= 0; i--)
	{
		LAPACKE_dormqr(LAPACK_COL_MAJOR, 'L', 'N', 2*n, N, n, &Q[(i + N*i)*n], N, tau + i*n, &H[i*n], N);
	}

	// entrywise absolute error (difference between Q*R and original H)
	double err = UniformDistance(N*N, H, H_ref);
	printf("Largest entrywise absolute error: %g\n", err);

	// clean up
	MKL_free(Q);
	MKL_free(H_ref);
	MKL_free(tau);
	MKL_free(H);

	return (err < 4e-16 ? 0 : 1);
}
