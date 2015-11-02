#include "measurement.h"
#include "random.h"
#include "util.h"
#include <math.h>
#include <stdio.h>


int MeasurementTest()
{
	int i, n;

	// lattice field dimensions
	#define Nx 4
	#define Ny 6
	// total number of lattice sites
	#define N  (Nx * Ny)

	// coupling constant in the Hubbard hamiltonian
	const double U = 6;

	// imaginary-time step size
	const double dt = 1.0/8;

	// t' (next-nearest neighbor) hopping parameter
	const double tp = -0.2;

	// chemical potential
	const double mu = -1.7;

	// number of time steps
	#define L 16

	// largest number of B_l matrices multiplied together before performing a QR decomposition
	const int prodBlen = 4;

	// number of measurement iterations
	const int nsampl = 12;

	printf("Comparing equal time measurements to unequal time measurements at time difference zero on a %i x %i lattice and beta = %g...\n", Nx, Ny, L*dt);

	// Hubbard-Stratonovich parameters
	stratonovich_params_t stratonovich_params;
	FillStratonovichParameters(U, dt, &stratonovich_params);
	
	// calculate matrix exponential of the kinetic nearest neighbor hopping matrix
	kinetic_t kinetic;
	SquareLatticeKineticExponential(Nx, Ny, tp, mu, dt, &kinetic);

	// artificial UNIX time
	time_t itime = 1420000241;

	// random generator seed; multiplicative constant from Pierre L'Ecuyer's paper
	randseed_t seed;
	Random_SeedInit(1865811235122147685LL * (uint64_t)itime, &seed);

	// time step matrices
	time_step_matrices_t tsm_u;
	time_step_matrices_t tsm_d;
	AllocateTimeStepMatrices(N, L, prodBlen, &tsm_u);
	AllocateTimeStepMatrices(N, L, prodBlen, &tsm_d);

	// allocate spin-up and spin-down Green's functions
	greens_func_t Gu, Gd;
	AllocateGreensFunction(N, &Gu);
	AllocateGreensFunction(N, &Gd);

	// allocate equal time measurement data structure
	measurement_data_t meas_data;
	AllocateMeasurementData(Nx, Ny, &meas_data);

	// allocate unequal time measurement data structure
	measurement_data_unequal_time_t meas_data_uneqlt;
	int status = AllocateUnequalTimeMeasurementData(Nx, Ny, L, &meas_data_uneqlt);
	if (status != 0) { return status; }

	// accumulate measurements of pseudo-random data
	for (n = 0; n < nsampl; n++)
	{
		// random Hubbard-Stratonovich field
		spin_field_t s[L*N];
		for (i = 0; i < L*N; i++)
		{
			s[i] = 	Random_GetUint(&seed) & 1;
		}

		InitTimeStepMatrices(&kinetic, stratonovich_params.expVu, s, &tsm_u);
		InitTimeStepMatrices(&kinetic, stratonovich_params.expVd, s, &tsm_d);

		GreenConstruct(&tsm_u, 0, &Gu);
		GreenConstruct(&tsm_d, 0, &Gd);

		AccumulateMeasurement(&Gu, &Gd, &meas_data);
		AccumulateUnequalTimeMeasurement((double)(Gu.sgndet * Gd.sgndet), tsm_u.B, tsm_d.B, &meas_data_uneqlt);
	}

	NormalizeMeasurementData(&meas_data);
	NormalizeUnequalTimeMeasurementData(&meas_data_uneqlt);

	// total density correlations (equal time)
	double nn_corr[N];
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
	err = fmax(err, fabs(meas_data.density_u - meas_data.uu_corr[0]));
	err = fmax(err, fabs(meas_data.density_d - meas_data.dd_corr[0]));
	err = fmax(err, fabs(meas_data.doubleocc - 0.5 * meas_data.ud_corr[0]));
	printf("Largest entrywise absolute error: %g\n", err);

	// clean up
	DeleteUnequalTimeMeasurementData(&meas_data_uneqlt);
	DeleteMeasurementData(&meas_data);
	DeleteGreensFunction(&Gd);
	DeleteGreensFunction(&Gu);
	DeleteTimeStepMatrices(&tsm_d);
	DeleteTimeStepMatrices(&tsm_u);
	DeleteKineticExponential(&kinetic);

	return (err < 5e-15 ? 0 : 1);
}
