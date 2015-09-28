#include "stratonovich.h"
#include "util.h"
#include <math.h>
#include <stdio.h>


int StratonovichTest()
{
	const double U = 4.0;
	const double dt = 1.0/11;

	// Hubbard-Stratonovich parameters
	printf("Calculating Hubbard-Stratonovich parameters for U = %g and dt = %g...\n", U, dt);
	stratonovich_params_t stratonovich_params;
	FillStratonovichParameters(U, dt, &stratonovich_params);

	// load reference data from disk
	stratonovich_params_t stratonovich_params_ref;
	ReadData("../test/stratonovich_test.dat", &stratonovich_params_ref, sizeof(stratonovich_params_ref), 1);

	// entrywise absolute error
	double err = 0;
	int i;
	for (i = 0; i < 2; i++)
	{
		err = fmax(err, fabs(stratonovich_params.expVu[i] - stratonovich_params_ref.expVu[i]));
		err = fmax(err, fabs(stratonovich_params.expVd[i] - stratonovich_params_ref.expVd[i]));
		err = fmax(err, fabs(stratonovich_params.delta[i] - stratonovich_params_ref.delta[i]));
	}
	printf("Largest absolute error: %g\n", err);

	return (err < 5e-16 ? 0 : 1);
}
