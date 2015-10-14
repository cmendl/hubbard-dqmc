#include "monte_carlo.h"
#include "util.h"
#include <mkl.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>


int MonteCarloPhononBlockTest()
{
	// lattice field dimensions
	#define Nx 4
	#define Ny 6
	// total number of lattice sites
	#define N  (Nx * Ny)

	// coupling constant in the Hubbard hamiltonian
	const double U = 4;

	// imaginary-time step size
	const double dt = 1.0/8;

	// number of time steps
	#define L 16

	// largest number of B_l matrices multiplied together before performing a QR decomposition
	const int prodBlen = 4;

	// Hubbard-Stratonovich parameters
	stratonovich_params_t stratonovich_params;
	FillStratonovichParameters(U, dt, &stratonovich_params);

	// phonon parameters
	phonon_params_t phonon_params;
	phonon_params.omega = 1.2;
	phonon_params.g = 0.65;
	phonon_params.box_width = 2;
	phonon_params.nblock_updates = 2;

	// calculate matrix exponential of the kinetic nearest neighbor hopping matrix
	kinetic_t kinetic;
	NearestNeighborKineticExponential(Nx, Ny, 0.0, dt, &kinetic);

	// Hubbard-Stratonovich field remains constant during phonon block updates
	spin_field_t s[L*N];
	int status;
	status = ReadData("../test/monte_carlo_phonon_block_test_s.dat", s, sizeof(spin_field_t), L*N); if (status != 0) { return status; }

	// initial phonon field
	double *X    = (double *)MKL_malloc(L*N * sizeof(double), MEM_DATA_ALIGN);
	double *expX = (double *)MKL_malloc(L*N * sizeof(double), MEM_DATA_ALIGN);
	status = ReadData("../test/monte_carlo_phonon_block_test_X0.dat", X, sizeof(double), L*N); if (status != 0) { return status; }
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

	// allocate and construct initial Green's functions
	greens_func_t Gu, Gd;
	AllocateGreensFunction(N, &Gu);
	AllocateGreensFunction(N, &Gd);
	GreenConstruct(&tsm_u, 0, &Gu);
	GreenConstruct(&tsm_d, 0, &Gd);

	// artificial UNIX time
	time_t itime = 1420000246;

	// random generator seed; multiplicative constant from Pierre L'Ecuyer's paper
	randseed_t seed;
	Random_SeedInit(1865811235122147685LL * (uint64_t)itime, &seed);

	// perform phonon block updates (first update will be accepted and second rejected)
	printf("Performing phonon block updates on a %i x %i lattice at beta = %g...\n", Nx, Ny, L*dt);
	PhononBlockUpdates(dt, &kinetic, &stratonovich_params, &phonon_params, &seed, s, X, expX, &tsm_u, &tsm_d, &Gu, &Gd);

	// load reference phonon field from disk
	double    X_ref[L*N];
	double expX_ref[L*N];
	status = ReadData("../test/monte_carlo_phonon_block_test_X1.dat",       X_ref, sizeof(double), L*N); if (status != 0) { return status; }
	status = ReadData("../test/monte_carlo_phonon_block_test_expX1.dat", expX_ref, sizeof(double), L*N); if (status != 0) { return status; }

	// entrywise absolute error of the phonon field
	double errX = fmax(
		UniformDistance(L*N,    X,    X_ref),
		UniformDistance(L*N, expX, expX_ref));
	printf("Largest entrywise absolute error of the phonon field: %g\n", errX);

	// load reference Green's function matrices from disk
	double Gu_mat_ref[N*N];
	double Gd_mat_ref[N*N];
	double detGu_ref, detGd_ref;
	status = ReadData("../test/monte_carlo_phonon_block_test_Gu1.dat",    Gu_mat_ref, sizeof(double), N*N); if (status != 0) { return status; }
	status = ReadData("../test/monte_carlo_phonon_block_test_Gd1.dat",    Gd_mat_ref, sizeof(double), N*N); if (status != 0) { return status; }
	status = ReadData("../test/monte_carlo_phonon_block_test_detGu1.dat", &detGu_ref, sizeof(double), 1);   if (status != 0) { return status; }
	status = ReadData("../test/monte_carlo_phonon_block_test_detGd1.dat", &detGd_ref, sizeof(double), 1);   if (status != 0) { return status; }

	// entrywise relative error of the Green's function matrices
	double errG_rel = 0;
	for (i = 0; i < N*N; i++)
	{
		errG_rel = fmax(errG_rel, fabs((Gu.mat[i] - Gu_mat_ref[i])/Gu_mat_ref[i]));
		errG_rel = fmax(errG_rel, fabs((Gd.mat[i] - Gd_mat_ref[i])/Gd_mat_ref[i]));
	}
	printf("Largest entrywise relative error of the Green's function matrices: %g\n", errG_rel);

	// entrywise absolute error of the Green's function matrices
	double errG_abs = fmax(
		UniformDistance(N*N, Gu.mat, Gu_mat_ref),
		UniformDistance(N*N, Gd.mat, Gd_mat_ref));
	printf("Largest entrywise absolute error of the Green's function matrices: %g\n", errG_abs);

	// relative error of determinant
	double detGu = Gu.sgndet * exp(Gu.logdet);
	double detGd = Gd.sgndet * exp(Gd.logdet);
	double err_detG = fmax(
		fabs((detGu - detGu_ref) / detGu_ref),
		fabs((detGd - detGd_ref) / detGd_ref));
	printf("Relative determinant error: %g\n", err_detG);

	// clean up
	DeleteGreensFunction(&Gd);
	DeleteGreensFunction(&Gu);
	DeleteTimeStepMatrices(&tsm_d);
	DeleteTimeStepMatrices(&tsm_u);
	MKL_free(expX);
	MKL_free(X);
	DeleteKineticExponential(&kinetic);

	return (errX < 5e-16 && errG_rel < 5e-11 && errG_abs < 2e-14 && err_detG < 8e-14 ? 0 : 1);
}
