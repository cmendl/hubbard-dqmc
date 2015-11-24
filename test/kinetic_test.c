#include "kinetic.h"
#include "util.h"
#include <mkl.h>
#include <math.h>
#include <stdio.h>


int KineticTest()
{
	// single-orbital simulation parameters
	sim_params_t params;
	AllocateSimulationParameters(1, &params);

	// lattice dimension
	params.Nx = 4;
	params.Ny = 6;
	// total number of lattice sites
	const int N = params.Nx * params.Ny;

	// imaginary-time step size
	params.dt = 1.0/7;

	// t (nearest neighbor) hopping parameter
	const double t = 1.0;
	params.t.ab[0] = t;
	params.t.ac[0] = t;

	// t' (next-nearest neighbor) hopping parameter
	const double tp = -2.0/13;
	params.t.ad[0] = tp;
	params.t.bc[0] = tp;

	// chemical potential
	params.mu = 2.0/9;

	// calculate matrix exponential of the kinetic nearest neighbor hopping matrix
	printf("Calculating matrix exponential of the kinetic nearest and next-nearest neighbor hopping matrix on a %i x %i lattice...\n", params.Nx, params.Ny);
	kinetic_t kinetic;
	RectangularKineticExponential(&params, &kinetic);

	// load reference data from disk
	double *expK_ref     = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	double *inv_expK_ref = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	ReadData("../test/kinetic_test_expK.dat",    expK_ref,     sizeof(double), N*N);
	ReadData("../test/kinetic_test_invexpK.dat", inv_expK_ref, sizeof(double), N*N);

	// compare with reference
	double err = fmax(
		UniformDistance(N*N, kinetic.expK,         expK_ref),
		UniformDistance(N*N, kinetic.inv_expK, inv_expK_ref));
	printf("Largest entrywise absolute error: %g\n", err);

	// clean up
	MKL_free(inv_expK_ref);
	MKL_free(expK_ref);
	DeleteKineticExponential(&kinetic);
	DeleteSimulationParameters(&params);

	return (err < 2e-15 ? 0 : 1);
}
