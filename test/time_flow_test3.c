#include "time_flow.h"
#include "kinetic.h"
#include "util.h"
#include <mkl.h>
#include <math.h>
#include <stdio.h>


int TimeFlowTest3()
{
	// lattice field dimensions
	#define Nx 4
	#define Ny 6
	// total number of lattice sites
	#define N  (Nx * Ny)

	// imaginary-time step size
	const double dt = 1.0/8;

	// electron-phonon interaction strength
	const double g = 0.7;

	// number of time steps
	#define L 16

	// largest number of B matrices multiplied together before performing a QR decomposition; must divide L
	const int prodBlen = 8;

	// lambda parameter (depends on U in Hamiltonian)
	const double lambda = 0.75;
	const double expV[2] = {
		exp(-lambda),
		exp( lambda)
	};

	// calculate matrix exponential of the kinetic nearest neighbor hopping matrix
	kinetic_t kinetic;
	NearestNeighborKineticExponential(Nx, Ny, 0.0, dt, &kinetic);

	// Hubbard-Stratonovich field
	const spin_field_t s[L*N] = {
		1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1,
		0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0,
		0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1,
		1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1,
		0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0,
		0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0,
		1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1,
		0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0,
		1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1,
		1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0,
		1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0,
		1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0,
		0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1,
		1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1,
		1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0,
		1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1
	};

	// phonon field
	double *X    = (double *)MKL_malloc(L*N * sizeof(double), MEM_DATA_ALIGN);
	double *expX = (double *)MKL_malloc(L*N * sizeof(double), MEM_DATA_ALIGN);
	ReadData("../test/time_flow_test3_X.dat", X, sizeof(double), L*N);
	int i;
	for (i = 0; i < L*N; i++)
	{
		expX[i] = exp(-dt*g * X[i]);
	}

	// allocate and initialize time step matrices
	time_step_matrices_t tsm;
	AllocateTimeStepMatrices(N, L, prodBlen, &tsm);
	InitPhononTimeStepMatrices(&kinetic, expV, s, expX, &tsm);

	// allocate structures for output, represented as Q * diag(d) * T
	double *Q   = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	double *tau = (double *)MKL_malloc(N   * sizeof(double), MEM_DATA_ALIGN);	// scalar factors of the elementary reflectors for the matrix Q
	double *d   = (double *)MKL_malloc(N   * sizeof(double), MEM_DATA_ALIGN);
	double *T   = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);

	// calculate the imaginary-time flow map
	printf("Calculating imaginary-time flow map including phonons on a %i x %i lattice...\n", Nx, Ny);
	TimeFlowMap(&tsm, 0, Q, tau, d, T);

	// form A = Q * diag(d) * T (only for comparison with reference)
	for (i = 0; i < N; i++)
	{
		cblas_dscal(N, d[i], &T[i], N);
	}
	// multiply by Q from the left
	LAPACKE_dormqr(LAPACK_COL_MAJOR, 'L', 'N', N, N, N, Q, N, tau, T, N);

	// load reference data from disk
	double A_ref[N*N];
	ReadData("../test/time_flow_test3_A.dat", A_ref, sizeof(double), N*N);

	// entrywise relative error
	double err_rel = 0;
	for (i = 0; i < N*N; i++)
	{
		err_rel = fmax(err_rel, fabs((T[i] - A_ref[i])/A_ref[i]));
	}
	printf("Largest entrywise relative error: %g\n", err_rel);

	// entrywise absolute error
	double err_abs = UniformDistance(N*N, T, A_ref);
	printf("Largest entrywise absolute error: %g\n", err_abs);

	// clean up
	MKL_free(T);
	MKL_free(d);
	MKL_free(tau);
	MKL_free(Q);
	DeleteTimeStepMatrices(&tsm);
	MKL_free(expX);
	MKL_free(X);
	DeleteKineticExponential(&kinetic);

	return (err_rel < 2e-14 && err_abs < 3e-9 ? 0 : 1);
}
