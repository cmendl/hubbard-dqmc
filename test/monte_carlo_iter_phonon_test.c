#include "monte_carlo.h"
#include "kinetic.h"
#include "stratonovich.h"
#include "time_flow.h"
#include "util.h"
#include <mkl.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>


void DQMCPhononIteration(const double dt, const kinetic_t *restrict kinetic, const stratonovich_params_t *restrict stratonovich_params, const phonon_params_t *restrict phonon_params, const int nwraps, randseed_t *restrict seed, spin_field_t *restrict s, double *restrict X, double *restrict expX, time_step_matrices_t *restrict tsm_u, time_step_matrices_t *restrict tsm_d, double *restrict Gu, double *restrict Gd);


int MonteCarloIterPhononTest()
{
	// lattice field dimensions
	#define Nx 4
	#define Ny 6
	// total number of lattice sites
	#define N  (Nx * Ny)

	// coupling constant in the Hubbard hamiltonian
	const double U = 4.5;

	// imaginary-time step size
	const double dt = 1.0/8;

	// number of time steps
	#define L 16

	// largest number of B_l matrices multiplied together before performing a QR decomposition
	const int prodBlen = 4;

	// number of "time slice wraps" before recomputing the Green's function
	const int nwraps = 8;

	// Hubbard-Stratonovich parameters
	stratonovich_params_t stratonovich_params;
	FillStratonovichParameters(U, dt, &stratonovich_params);

	// phonon parameters
	phonon_params_t phonon_params;
	phonon_params.omega = 1.3;
	phonon_params.g = 0.7;
	phonon_params.box_width = 12;
	phonon_params.nblock_updates = 0;	// disable block updates

	// calculate matrix exponential of the kinetic nearest neighbor hopping matrix
	kinetic_t kinetic;
	NearestNeighborKineticExponential(Nx, Ny, 0.0, dt, &kinetic);

	// initial Hubbard-Stratonovich field
	spin_field_t s[L*N] = {
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

	// initial phonon field
	double *X    = (double *)MKL_malloc(L*N * sizeof(double), MEM_DATA_ALIGN);
	double *expX = (double *)MKL_malloc(L*N * sizeof(double), MEM_DATA_ALIGN);
	int status;
	status = ReadData("../test/monte_carlo_iter_phonon_test_X0.dat", X, sizeof(double), L*N); if (status != 0) { return status; }
	int i;
	for (i = 0; i < L*N; i++)
	{
		expX[i] = exp(-dt*phonon_params.g * X[i]);
	}

	// time step matrices
	time_step_matrices_t tsm_u;
	time_step_matrices_t tsm_d;
	AllocateTimeStepMatrices(N, L, prodBlen, &tsm_u);
	AllocateTimeStepMatrices(N, L, prodBlen, &tsm_d);
	InitPhononTimeStepMatrices(&kinetic, stratonovich_params.expVu, s, expX, &tsm_u);
	InitPhononTimeStepMatrices(&kinetic, stratonovich_params.expVd, s, expX, &tsm_d);

	// spin-up and spin-down Green's function matrices; will be initialized during call of 'DQMCPhononIteration'
	double *Gu = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	double *Gd = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);

	// artificial UNIX time
	time_t itime = 1420000241;

	// random generator seed; multiplicative constant from Pierre L'Ecuyer's paper
	randseed_t seed;
	Random_SeedInit(1865811235122147685LL * (uint64_t)itime, &seed);

	// perform a Determinant Quantum Monte Carlo (DQMC) iteration
	printf("Performing a Determinant Quantum Monte Carlo (DQMC) iteration on a %i x %i lattice at beta = %g, taking phonons into account...\n", Nx, Ny, L*dt);
	DQMCPhononIteration(dt, &kinetic, &stratonovich_params, &phonon_params, nwraps, &seed, s, X, expX, &tsm_u, &tsm_d, Gu, Gd);

	// reference Hubbard-Stratonovich field after DQMC iteration
	spin_field_t s_ref[L*N] = {
		0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0,
		1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1,
		1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1,
		0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0,
		1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1,
		1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1,
		0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0,
		1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1,
		1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0,
		1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1,
		0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1,
		1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1,
		1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1,
		1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1,
		1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1,
		1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0
	};

	// number of deviating Hubbard-Stratonovich field entries
	int err_field = 0;
	for (i = 0; i < L*N; i++)
	{
		err_field += abs(s[i] - s_ref[i]);
	}
	printf("Number of deviating Hubbard-Stratonovich field entries: %i\n", err_field);

	// load reference phonon field from disk
	double X_ref[L*N];
	status = ReadData("../test/monte_carlo_iter_phonon_test_X1.dat", X_ref, sizeof(double), L*N); if (status != 0) { return status; }

	// entrywise absolute error of the phonon field
	double errX = 0;
	for (i = 0; i < L*N; i++)
	{
		errX = fmax(errX, fabs(X[i] - X_ref[i]));
	}
	printf("Largest entrywise absolute error of the phonon field: %g\n", errX);

	// load reference Green's functions from disk
	double Gu_ref[N*N];
	double Gd_ref[N*N];
	status = ReadData("../test/monte_carlo_iter_phonon_test_Gu1.dat", Gu_ref, sizeof(double), N*N); if (status != 0) { return status; }
	status = ReadData("../test/monte_carlo_iter_phonon_test_Gd1.dat", Gd_ref, sizeof(double), N*N); if (status != 0) { return status; }

	// entrywise relative error of the Green's function matrices
	double errG_rel = 0;
	for (i = 0; i < N*N; i++)
	{
		errG_rel = fmax(errG_rel, fabs((Gu[i] - Gu_ref[i])/Gu_ref[i]));
		errG_rel = fmax(errG_rel, fabs((Gd[i] - Gd_ref[i])/Gd_ref[i]));
	}
	printf("Largest entrywise relative error of the Green's function matrices: %g\n", errG_rel);

	// entrywise absolute error of the Green's function matrices
	double errG_abs = fmax(
		UniformDistance(N*N, Gu, Gu_ref),
		UniformDistance(N*N, Gd, Gd_ref));
	printf("Largest entrywise absolute error of the Green's function matrices: %g\n", errG_abs);

	// clean up
	MKL_free(Gd);
	MKL_free(Gu);
	DeleteTimeStepMatrices(&tsm_d);
	DeleteTimeStepMatrices(&tsm_u);
	MKL_free(expX);
	MKL_free(X);
	DeleteKineticExponential(&kinetic);

	return (err_field == 0 && errX < 2e-15 && errG_rel < 4e-8 && errG_abs < 4e-11 ? 0 : 1);
}
