#include "kinetic.h"
#include "util.h"
#include <mkl.h>
#include <math.h>
#include <stdio.h>


int KineticTest4()
{
	// three-orbital simulation parameters
	sim_params_t params = { 0 };
	AllocateSimulationParameters(3, &params);

	// lattice dimension
	params.Nx = 5;
	params.Ny = 6;
	// shift of x-coordinate when wrapping around in y-direction
	params.pbc_shift = params.Ny / 2;

	// t (nearest neighbor) hopping parameter
	const double t = 1.2;
	// set hoppings corresponding to a Kagome lattice
	params.t.aa[1] = t;		// (1,0) entry
	params.t.aa[2] = t;		// (2,0) entry
	params.t.aa[7] = t;		// (1,2) entry
	params.t.ac[1] = t;		// (1,0) entry
	params.t.bc[2] = t;		// (2,0) entry
	params.t.ab[7] = t;		// (1,2) entry

	// chemical potential
	params.mu = 2.0/7;

	// imaginary-time step size
	params.dt = 1.0/9;

	// calculate matrix exponential of the kinetic nearest neighbor hopping matrix
	printf("Calculating matrix exponential of the kinetic hopping matrix on a Kagome lattice...\n");
	kinetic_t kinetic;
	RectangularKineticExponential(&params, &kinetic);

	// load reference data from disk
	const int N = kinetic.Ncell * kinetic.Norb;
	double *expK_ref     = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	double *inv_expK_ref = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	ReadData("../test/kinetic_test4_expK.dat",    expK_ref,     sizeof(double), N*N);
	ReadData("../test/kinetic_test4_invexpK.dat", inv_expK_ref, sizeof(double), N*N);

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
