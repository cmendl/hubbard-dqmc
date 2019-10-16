#include "measurement.h"
#include "random.h"
#include "profiler.h"
#include "util.h"
#include <math.h>
#include <stdio.h>


int MeasurementTest()
{
	int i, n;
	sim_params_t params = { 0 };
	AllocateSimulationParameters(1, &params);

	// lattice field dimensions
	params.Nx = 4;
	params.Ny = 6;

	// without shift when wrapping around in y-direction
	params.pbc_shift = 0;

	// total number of lattice sites
	const int N = params.Norb * params.Nx*params.Ny;

	// coupling constant in the Hubbard hamiltonian
	params.U[0] = 6.0;

	// imaginary-time step size
	params.dt = 1.0/8;

	// hopping parameters
	params.t.ab[0] = 1;
	params.t.ac[0] = 1;
	params.t.ad[0] = -0.2;
	params.t.bc[0] = -0.2;

	// chemical potential
	params.mu = -1.7;
	params.eps[0] = 0;

	// number of time steps
	params.L = 16;

	// largest number of B_l matrices multiplied together before performing a QR decomposition
	params.prodBlen = 4;

	// number of measurement iterations
	params.nsampl = 12;

	// initialize profiling (to avoid runtime exception: profiler called during GreenConstruct)
	Profile_Start();

	printf("Comparing equal time measurements to unequal time measurements at time difference zero on a %i x %i lattice and beta = %g...\n", params.Nx, params.Ny, params.L*params.dt);

	// Hubbard-Stratonovich parameters
	stratonovich_params_t stratonovich_params;
	FillStratonovichParameters(params.Norb, params.U, params.dt, &stratonovich_params);

	// calculate matrix exponential of the kinetic nearest neighbor hopping matrix
	kinetic_t kinetic;
	RectangularKineticExponential(&params, &kinetic);

	// artificial UNIX time
	params.itime = 1420000241;

	// random generator seed; multiplicative constant from Pierre L'Ecuyer's paper
	randseed_t seed;
	Random_SeedInit(1865811235122147685LL * params.itime, &seed);

	// time step matrices
	time_step_matrices_t tsm_u;
	time_step_matrices_t tsm_d;
	AllocateTimeStepMatrices(N, params.L, params.prodBlen, &tsm_u);
	AllocateTimeStepMatrices(N, params.L, params.prodBlen, &tsm_d);

	// allocate spin-up and spin-down Green's functions
	greens_func_t Gu, Gd;
	AllocateGreensFunction(N, &Gu);
	AllocateGreensFunction(N, &Gd);

	// allocate equal time measurement data structure
	measurement_data_t meas_data;
	AllocateMeasurementData(params.Norb, params.Nx, params.Ny, params.pbc_shift, &meas_data);

	// allocate unequal time measurement data structure
	measurement_data_unequal_time_t meas_data_uneqlt;
	int status = AllocateUnequalTimeMeasurementData(params.Norb, params.Nx, params.Ny, params.pbc_shift, params.L, tsm_u.numBprod, &meas_data_uneqlt);
	if (status != 0) { return status; }

	// accumulate measurements of pseudo-random data
	spin_field_t *s = (spin_field_t *)algn_malloc(params.L*N * sizeof(spin_field_t));
	for (n = 0; n < params.nsampl; n++)
	{
		// random Hubbard-Stratonovich field
		for (i = 0; i < params.L*N; i++)
		{
			s[i] = 	Random_GetUint(&seed) & 1;
		}

		InitTimeStepMatrices(&kinetic, stratonovich_params.expVu, s, &tsm_u);
		InitTimeStepMatrices(&kinetic, stratonovich_params.expVd, s, &tsm_d);

		GreenConstruct(&tsm_u, 0, &Gu);
		GreenConstruct(&tsm_d, 0, &Gd);

		AccumulateMeasurement(&Gu, &Gd, &meas_data);
		AccumulateUnequalTimeMeasurement((double)(Gu.sgndet * Gd.sgndet), &tsm_u, &tsm_d, &meas_data_uneqlt);
	}

	NormalizeMeasurementData(&meas_data);
	NormalizeUnequalTimeMeasurementData(&meas_data_uneqlt);

	// total density correlations (equal time)
	double *nn_corr = (double *)algn_malloc(N * sizeof(double));
	for (i = 0; i < N; i++)
	{
		nn_corr[i] = meas_data.uu_corr[i] + meas_data.dd_corr[i] + meas_data.ud_corr[i];
	}

	// compare
	double err = 0;
	err = fmax(err, UniformDistance(N, meas_data_uneqlt.nn_corr,           nn_corr));
	err = fmax(err, UniformDistance(N, meas_data_uneqlt.zz_corr, meas_data.zz_corr));
	err = fmax(err, UniformDistance(N, meas_data_uneqlt.xx_corr, meas_data.xx_corr));
	err = fmax(err, fabs(meas_data_uneqlt.sign - meas_data.sign));
	err = fmax(err, fabs(meas_data.density_u[0] - meas_data.uu_corr[0]));
	err = fmax(err, fabs(meas_data.density_d[0] - meas_data.dd_corr[0]));
	err = fmax(err, fabs(meas_data.doubleocc[0] - 0.5 * meas_data.ud_corr[0]));
	printf("Largest entrywise absolute error: %g\n", err);

	// clean up
	Profile_Stop();
	algn_free(nn_corr);
	algn_free(s);
	DeleteUnequalTimeMeasurementData(&meas_data_uneqlt);
	DeleteMeasurementData(&meas_data);
	DeleteGreensFunction(&Gd);
	DeleteGreensFunction(&Gu);
	DeleteTimeStepMatrices(&tsm_d);
	DeleteTimeStepMatrices(&tsm_u);
	DeleteKineticExponential(&kinetic);
	DeleteStratonovichParameters(&stratonovich_params);
	DeleteSimulationParameters(&params);

	return (err < 5e-15 ? 0 : 1);
}
