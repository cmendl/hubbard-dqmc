#include "monte_carlo.h"
#include "util.h"
#include <mkl.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>


int MonteCarloIterTest()
{
	// two-orbital simulation parameters
	sim_params_t params = { 0 };
	AllocateSimulationParameters(2, &params);

	// lattice field dimensions
	params.Nx = 4;
	params.Ny = 5;

	// total number of lattice sites
	const int N = params.Norb * params.Nx*params.Ny;

	// coupling constants in the Hubbard hamiltonian
	params.U[0] = 4.5;
	params.U[1] = 11.0/3;

	// imaginary-time step size
	params.dt = 1.0/8;

	int status;

	// read hopping parameters from disk
	status = ReadData("../test/monte_carlo_iter_test_taa.dat", params.t.aa, sizeof(double), params.Norb*params.Norb);	if (status < 0) { return status; }
	status = ReadData("../test/monte_carlo_iter_test_tab.dat", params.t.ab, sizeof(double), params.Norb*params.Norb);	if (status < 0) { return status; }
	status = ReadData("../test/monte_carlo_iter_test_tac.dat", params.t.ac, sizeof(double), params.Norb*params.Norb);	if (status < 0) { return status; }
	status = ReadData("../test/monte_carlo_iter_test_tad.dat", params.t.ad, sizeof(double), params.Norb*params.Norb);	if (status < 0) { return status; }
	status = ReadData("../test/monte_carlo_iter_test_tbc.dat", params.t.bc, sizeof(double), params.Norb*params.Norb);	if (status < 0) { return status; }

	// chemical potential
	params.mu = 5.0/7;

	// site energies
	params.eps[0] = -1.0/5;
	params.eps[1] =  2.0/11;

	// number of time steps
	params.L = 16;

	// largest number of B_l matrices multiplied together before performing a QR decomposition
	params.prodBlen = 4;

	// number of "time slice wraps" before recomputing the Green's function
	params.nwraps = 8;

	// Hubbard-Stratonovich parameters
	stratonovich_params_t stratonovich_params;
	FillStratonovichParameters(params.Norb, params.U, params.dt, &stratonovich_params);

	// calculate matrix exponential of the kinetic nearest neighbor hopping matrix
	kinetic_t kinetic;
	RectangularKineticExponential(&params, &kinetic);

	// load initial Hubbard-Stratonovich field from disk
	spin_field_t *s = (spin_field_t *)MKL_malloc(N*params.L *sizeof(spin_field_t), MEM_DATA_ALIGN);
	status = ReadData("../test/monte_carlo_iter_test_HS0.dat", s, sizeof(spin_field_t), N*params.L);
	if (status != 0) { return status; }

	// time step matrices
	time_step_matrices_t tsm_u;
	time_step_matrices_t tsm_d;
	AllocateTimeStepMatrices(N, params.L, params.prodBlen, &tsm_u);
	AllocateTimeStepMatrices(N, params.L, params.prodBlen, &tsm_d);
	InitTimeStepMatrices(&kinetic, stratonovich_params.expVu, s, &tsm_u);
	InitTimeStepMatrices(&kinetic, stratonovich_params.expVd, s, &tsm_d);

	// allocate and construct spin-up and spin-down Green's functions
	greens_func_t Gu, Gd;
	AllocateGreensFunction(N, &Gu);
	AllocateGreensFunction(N, &Gd);
	GreenConstruct(&tsm_u, 0, &Gu);
	GreenConstruct(&tsm_d, 0, &Gd);

	// artificial UNIX time
	params.itime = 1420000241;

	// random generator seed; multiplicative constant from Pierre L'Ecuyer's paper
	randseed_t seed;
	Random_SeedInit(1865811235122147685LL * params.itime, &seed);

	// perform a Determinant Quantum Monte Carlo (DQMC) iteration
	printf("Performing a Determinant Quantum Monte Carlo (DQMC) iteration on a %i x %i lattice with %i orbitals per unit cell at beta = %g...\n", params.Nx, params.Ny, params.Norb, params.L*params.dt);
	DQMCIteration(&kinetic, &stratonovich_params, params.nwraps, &seed, s, &tsm_u, &tsm_d, &Gu, &Gd);

	// reference Hubbard-Stratonovich field after DQMC iteration
	spin_field_t *s_ref = (spin_field_t *)MKL_malloc(N*params.L *sizeof(spin_field_t), MEM_DATA_ALIGN);
	status = ReadData("../test/monte_carlo_iter_test_HS1.dat", s_ref, sizeof(spin_field_t), N*params.L);
	if (status != 0) { return status; }

	// number of deviating Hubbard-Stratonovich field entries
	int err_field = 0;
	int i;
	for (i = 0; i < params.L*N; i++)
	{
		err_field += abs(s[i] - s_ref[i]);
	}
	printf("Number of deviating Hubbard-Stratonovich field entries: %i\n", err_field);

	// load reference data from disk
	double *Gu_mat_ref = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	double *Gd_mat_ref = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	double detGu_ref, detGd_ref;
	status = ReadData("../test/monte_carlo_iter_test_Gu1.dat",    Gu_mat_ref, sizeof(double), N*N); if (status != 0) { return status; }
	status = ReadData("../test/monte_carlo_iter_test_Gd1.dat",    Gd_mat_ref, sizeof(double), N*N); if (status != 0) { return status; }
	status = ReadData("../test/monte_carlo_iter_test_detGu1.dat", &detGu_ref, sizeof(double), 1);   if (status != 0) { return status; }
	status = ReadData("../test/monte_carlo_iter_test_detGd1.dat", &detGd_ref, sizeof(double), 1);   if (status != 0) { return status; }

	// entrywise relative error of the Green's function matrices
	double err_rel = 0;
	for (i = 0; i < N*N; i++)
	{
		err_rel = fmax(err_rel, fabs((Gu.mat[i] - Gu_mat_ref[i])/Gu_mat_ref[i]));
		err_rel = fmax(err_rel, fabs((Gd.mat[i] - Gd_mat_ref[i])/Gd_mat_ref[i]));
	}
	printf("Largest entrywise relative error of the Green's function matrices: %g\n", err_rel);

	// entrywise absolute error of the Green's function matrices
	double err_abs = fmax(
		UniformDistance(N*N, Gu.mat, Gu_mat_ref),
		UniformDistance(N*N, Gd.mat, Gd_mat_ref));
	printf("Largest entrywise absolute error of the Green's function matrices: %g\n", err_abs);

	// relative error of determinant
	double detGu = Gu.sgndet * exp(Gu.logdet);
	double detGd = Gd.sgndet * exp(Gd.logdet);
	double err_det = fmax(
		fabs((detGu - detGu_ref) / detGu_ref),
		fabs((detGd - detGd_ref) / detGd_ref));
	printf("Relative determinant error: %g\n", err_det);

	// clean up
	MKL_free(Gd_mat_ref);
	MKL_free(Gu_mat_ref);
	MKL_free(s_ref);
	DeleteGreensFunction(&Gd);
	DeleteGreensFunction(&Gu);
	DeleteTimeStepMatrices(&tsm_d);
	DeleteTimeStepMatrices(&tsm_u);
	MKL_free(s);
	DeleteKineticExponential(&kinetic);
	DeleteStratonovichParameters(&stratonovich_params);
	DeleteSimulationParameters(&params);

	return (err_field == 0 && err_rel < 1e-10 && err_abs < 4e-14 && err_det < 2e-13 ? 0 : 1);
}
