#include "time_flow.h"
#include "kinetic.h"
#include "profiler.h"
#include "util.h"
#include <mkl.h>
#include <math.h>
#include <stdio.h>


int TimeFlowTest2()
{
	// single-orbital simulation parameters
	sim_params_t params = { 0 };
	AllocateSimulationParameters(1, &params);

	// lattice field dimensions
	params.Nx = 5;
	params.Ny = 6;
	// total number of lattice sites
	const int N = params.Nx * params.Ny;

	// imaginary-time step size
	params.dt = 1.0/11;

	// t (nearest neighbor) hopping parameter
	const double t = 1.0;
	params.t.ab[0] = t;
	params.t.ac[0] = t;

	// chemical potential
	params.mu = 0;

	// number of time steps
	params.L = 21;

	// largest number of B matrices multiplied together before performing a QR decomposition; must divide L
	params.prodBlen = 7;

	// lambda parameter (depends on U in Hamiltonian)
	const double lambda = 0.75;

	// calculate matrix exponential of the kinetic nearest neighbor hopping matrix
	kinetic_t kinetic;
	RectangularKineticExponential(&params, &kinetic);

	// Hubbard-Stratonovich field
	const spin_field_t s[] = {
		1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1,
		0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1,
		0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0,
		1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0,
		0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1,
		0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1,
		1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0,
		0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0,
		1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1,
		1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0,
		1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0,
		1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1,
		0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0,
		1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0,
		1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0,
		1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0,
		1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1,
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0,
		0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0,
		1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1,
		1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0
	};

	double *expV[2] = { (double *)MKL_malloc(sizeof(double), MEM_DATA_ALIGN), (double *)MKL_malloc(sizeof(double), MEM_DATA_ALIGN) };
	expV[0][0] = exp(-lambda);
	expV[1][0] = exp( lambda);

	// allocate and initialize time step matrices
	time_step_matrices_t tsm;
	AllocateTimeStepMatrices(N, params.L, params.prodBlen, &tsm);
	InitTimeStepMatrices(&kinetic, expV, s, &tsm);

	// allocate structures for output, represented as Q * diag(d) * T
	double *Q   = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	double *tau = (double *)MKL_malloc(N   * sizeof(double), MEM_DATA_ALIGN);	// scalar factors of the elementary reflectors for the matrix Q
	double *d   = (double *)MKL_malloc(N   * sizeof(double), MEM_DATA_ALIGN);
	double *T   = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);

	// initialize profiling (to avoid runtime exception: profiler called by TimeFlowMap)
	Profile_Start();

	// calculate the imaginary-time flow map
	printf("Calculating imaginary-time flow map on a %i x %i lattice with a single orbital per unit cell...\n", params.Nx, params.Ny);
	TimeFlowMap(&tsm, 0, Q, tau, d, T);

	// form A = Q * diag(d) * T (only for comparison with reference)
	int i;
	for (i = 0; i < N; i++)
	{
		cblas_dscal(N, d[i], &T[i], N);
	}
	// multiply by Q from the left
	LAPACKE_dormqr(LAPACK_COL_MAJOR, 'L', 'N', N, N, N, Q, N, tau, T, N);

	// load reference data from disk
	double *A_ref = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	ReadData("../test/time_flow_test2.dat", A_ref, sizeof(double), N*N);

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
	Profile_Stop();
	MKL_free(A_ref);
	MKL_free(T);
	MKL_free(d);
	MKL_free(tau);
	MKL_free(Q);
	DeleteTimeStepMatrices(&tsm);
	MKL_free(expV[1]);
	MKL_free(expV[0]);
	DeleteKineticExponential(&kinetic);
	DeleteSimulationParameters(&params);

	return (err_rel < 2e-13 && err_abs < 4e-8 ? 0 : 1);
}
