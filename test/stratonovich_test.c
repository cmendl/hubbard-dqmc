#include "stratonovich.h"
#include "util.h"
#include <mkl.h>
#include <math.h>
#include <stdio.h>


int StratonovichTest()
{
	const int Norb = 4;
	const double U[] = { 4.2, 2.7, 5.1, 3.8 };
	const double dt = 1.0/11;

	// Hubbard-Stratonovich parameters
	printf("Calculating Hubbard-Stratonovich parameters for %i orbitals and dt = %g...\n", Norb, dt);
	stratonovich_params_t stratonovich_params;
	FillStratonovichParameters(Norb, U, dt, &stratonovich_params);

	// load reference data from disk
	double *expVu_ref = (double *)MKL_malloc(2*Norb * sizeof(double), MEM_DATA_ALIGN);
	double *expVd_ref = (double *)MKL_malloc(2*Norb * sizeof(double), MEM_DATA_ALIGN);
	double *delta_ref = (double *)MKL_malloc(2*Norb * sizeof(double), MEM_DATA_ALIGN);
	ReadData("../test/stratonovich_test_expVu.dat", expVu_ref, sizeof(double), 2*Norb);
	ReadData("../test/stratonovich_test_expVd.dat", expVd_ref, sizeof(double), 2*Norb);
	ReadData("../test/stratonovich_test_delta.dat", delta_ref, sizeof(double), 2*Norb);

	// entrywise absolute error
	double err = 0;
	int i;
	for (i = 0; i < 2; i++)
	{
		int o;
		for (o = 0; o < Norb; o++)
		{

			err = fmax(err, fabs(stratonovich_params.expVu[i][o] - expVu_ref[o + i*Norb]));
			err = fmax(err, fabs(stratonovich_params.expVd[i][o] - expVd_ref[o + i*Norb]));
			err = fmax(err, fabs(stratonovich_params.delta[i][o] - delta_ref[o + i*Norb]));
		}
	}
	printf("Largest absolute error: %g\n", err);

	// clean up
	MKL_free(delta_ref);
	MKL_free(expVd_ref);
	MKL_free(expVu_ref);
	DeleteStratonovichParameters(&stratonovich_params);

	return (err < 5e-16 ? 0 : 1);
}
