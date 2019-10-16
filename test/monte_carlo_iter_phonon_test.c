#include "monte_carlo.h"
#include "profiler.h"
#include "util.h"
#include <math.h>
#include <assert.h>
#include <stdio.h>


int MonteCarloIterPhononTest()
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
	params.U[0] = 10.0/3;
	params.U[1] = 11.0/5;

	// imaginary-time step size
	params.dt = 1.0/8;

	int status;

	// read hopping parameters from disk
	status = ReadData("../test/monte_carlo_iter_phonon_test_taa.dat", params.t.aa, sizeof(double), params.Norb*params.Norb);  if (status < 0) { return status; }
	status = ReadData("../test/monte_carlo_iter_phonon_test_tab.dat", params.t.ab, sizeof(double), params.Norb*params.Norb);  if (status < 0) { return status; }
	status = ReadData("../test/monte_carlo_iter_phonon_test_tac.dat", params.t.ac, sizeof(double), params.Norb*params.Norb);  if (status < 0) { return status; }
	status = ReadData("../test/monte_carlo_iter_phonon_test_tad.dat", params.t.ad, sizeof(double), params.Norb*params.Norb);  if (status < 0) { return status; }
	status = ReadData("../test/monte_carlo_iter_phonon_test_tbc.dat", params.t.bc, sizeof(double), params.Norb*params.Norb);  if (status < 0) { return status; }

	// chemical potential
	params.mu = 6.0/7;

	// site energies
	params.eps[0] = -1.0/5;
	params.eps[1] =  3.0/11;

	// number of time steps
	params.L = 16;

	// largest number of B_l matrices multiplied together before performing a QR decomposition
	params.prodBlen = 4;

	// number of "time slice wraps" before recomputing the Green's function
	params.nwraps = 8;

	// Hubbard-Stratonovich parameters
	stratonovich_params_t stratonovich_params;
	FillStratonovichParameters(params.Norb, params.U, params.dt, &stratonovich_params);

	// phonon parameters
	params.phonon_params.omega[0] = 1.3;
	params.phonon_params.omega[1] = 7.0/8;
	params.phonon_params.g[0] = 5.0/11;
	params.phonon_params.g[1] = 7.0/10;
	params.phonon_params.local_box_width = 12;
	params.phonon_params.n_local_updates = 40;
	params.phonon_params.n_block_updates = 0;   // disable block updates


	// initialize profiling (to avoid runtime exception: profiler called by DQMCPhononIteration)
	Profile_Start();

	// calculate matrix exponential of the kinetic nearest neighbor hopping matrix
	kinetic_t kinetic;
	RectangularKineticExponential(&params, &kinetic);

	// load initial Hubbard-Stratonovich field from disk
	spin_field_t *s = (spin_field_t *)algn_malloc(N*params.L *sizeof(spin_field_t));
	status = ReadData("../test/monte_carlo_iter_phonon_test_HS0.dat", s, sizeof(spin_field_t), N*params.L);  if (status != 0) { return status; }

	// initial phonon field
	double *X    = (double *)algn_malloc(params.L*N * sizeof(double));
	double *expX = (double *)algn_malloc(params.L*N * sizeof(double));
	status = ReadData("../test/monte_carlo_iter_phonon_test_X0.dat", X, sizeof(double), params.L*N);  if (status != 0) { return status; }
	int i, l;
	for (l = 0; l < params.L; l++)
	{
		for (i = 0; i < N; i++)
		{
			const int o = i / kinetic.Ncell;    // orbital index
			expX[i + l*N] = exp(-params.dt*params.phonon_params.g[o] * X[i + l*N]);
		}
	}

	// time step matrices
	time_step_matrices_t tsm_u;
	time_step_matrices_t tsm_d;
	AllocateTimeStepMatrices(N, params.L, params.prodBlen, &tsm_u);
	AllocateTimeStepMatrices(N, params.L, params.prodBlen, &tsm_d);
	InitPhononTimeStepMatrices(&kinetic, stratonovich_params.expVu, s, expX, &tsm_u);
	InitPhononTimeStepMatrices(&kinetic, stratonovich_params.expVd, s, expX, &tsm_d);

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

	// phonon measurement data (unused here)
	measurement_data_phonon_t meas_data_phonon;
	AllocatePhononData(params.Norb, params.Nx, params.Ny, params.pbc_shift, params.L, 1, &meas_data_phonon);

	// perform a Determinant Quantum Monte Carlo (DQMC) iteration
	printf("Performing a Determinant Quantum Monte Carlo (DQMC) iteration on a %i x %i lattice with %i orbitals per unit cell, taking phonons into account...\n", params.Nx, params.Ny, params.Norb);
	DQMCPhononIteration(params.dt, params.mu, &kinetic, false, &stratonovich_params, &params.phonon_params, params.nwraps, &seed, s, X, expX, &tsm_u, &tsm_d, &Gu, &Gd, 0, NULL, &meas_data_phonon);

	// reference Hubbard-Stratonovich field after DQMC iteration
	spin_field_t *s_ref = (spin_field_t *)algn_malloc(N*params.L *sizeof(spin_field_t));
	status = ReadData("../test/monte_carlo_iter_phonon_test_HS1.dat", s_ref, sizeof(spin_field_t), N*params.L);  if (status != 0) { return status; }

	// number of deviating Hubbard-Stratonovich field entries
	int err_field = 0;
	for (i = 0; i < params.L*N; i++)
	{
		err_field += abs(s[i] - s_ref[i]);
	}
	printf("Number of deviating Hubbard-Stratonovich field entries: %i\n", err_field);

	// load reference phonon field from disk
	double *X_ref = (double *)algn_malloc(params.L*N * sizeof(double));
	status = ReadData("../test/monte_carlo_iter_phonon_test_X1.dat", X_ref, sizeof(double), params.L*N);  if (status != 0) { return status; }

	// entrywise absolute error of the phonon field
	double errX = 0;
	for (i = 0; i < params.L*N; i++)
	{
		errX = fmax(errX, fabs(X[i] - X_ref[i]));
	}
	printf("Largest entrywise absolute error of the phonon field: %g\n", errX);

	// load reference Green's functions from disk
	double *Gu_mat_ref = (double *)algn_malloc(N*N * sizeof(double));
	double *Gd_mat_ref = (double *)algn_malloc(N*N * sizeof(double));
	double detGu_ref, detGd_ref;
	status = ReadData("../test/monte_carlo_iter_phonon_test_Gu1.dat",    Gu_mat_ref, sizeof(double), N*N);  if (status != 0) { return status; }
	status = ReadData("../test/monte_carlo_iter_phonon_test_Gd1.dat",    Gd_mat_ref, sizeof(double), N*N);  if (status != 0) { return status; }
	status = ReadData("../test/monte_carlo_iter_phonon_test_detGu1.dat", &detGu_ref, sizeof(double), 1);    if (status != 0) { return status; }
	status = ReadData("../test/monte_carlo_iter_phonon_test_detGd1.dat", &detGd_ref, sizeof(double), 1);    if (status != 0) { return status; }

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
	Profile_Stop();
	algn_free(Gd_mat_ref);
	algn_free(Gu_mat_ref);
	algn_free(X_ref);
	algn_free(s_ref);
	DeletePhononData(&meas_data_phonon);
	DeleteGreensFunction(&Gd);
	DeleteGreensFunction(&Gu);
	DeleteTimeStepMatrices(&tsm_d);
	DeleteTimeStepMatrices(&tsm_u);
	algn_free(expX);
	algn_free(X);
	algn_free(s);
	DeleteKineticExponential(&kinetic);
	DeleteStratonovichParameters(&stratonovich_params);
	DeleteSimulationParameters(&params);

	return (err_field == 0 && errX < 1e-15 && errG_rel < 5e-11 && errG_abs < 1e-14 && err_detG < 4e-13 ? 0 : 1);
}
