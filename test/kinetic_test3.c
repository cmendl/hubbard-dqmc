#include "kinetic.h"
#include "util.h"
#include <math.h>
#include <stdio.h>


int KineticTest3()
{
	// two-orbital simulation parameters
	sim_params_t params = { 0 };
	AllocateSimulationParameters(2, &params);

	// lattice dimension
	params.Nx = 5;
	params.Ny = 6;
	// shift of x-coordinate when wrapping around in y-direction
	params.pbc_shift = params.Ny / 2;

	// t (nearest neighbor) hopping parameter
	const double t = 0.8;
	params.t.aa[1] = t;		// (1,0) entry
	params.t.ab[1] = t;		// (1,0) entry
	params.t.ac[1] = t;		// (1,0) entry

	// chemical potential
	params.mu = 2.0/9;

	// imaginary-time step size
	params.dt = 1.0/7;

	// calculate matrix exponential of the kinetic nearest neighbor hopping matrix
	printf("Calculating matrix exponential of the kinetic hopping matrix on a hexagonal lattice...\n");
	kinetic_t kinetic;
	RectangularKineticExponential(&params, &kinetic);

	// load reference data from disk
	const int N = kinetic.Ncell * kinetic.Norb;
	double *expK_ref     = (double *)algn_malloc(N*N * sizeof(double));
	double *inv_expK_ref = (double *)algn_malloc(N*N * sizeof(double));
	ReadData("../test/kinetic_test3_expK.dat",    expK_ref,     sizeof(double), N*N);
	ReadData("../test/kinetic_test3_invexpK.dat", inv_expK_ref, sizeof(double), N*N);

	// compare with reference
	double err = fmax(
		UniformDistance(N*N, kinetic.expK,         expK_ref),
		UniformDistance(N*N, kinetic.inv_expK, inv_expK_ref));
	printf("Largest entrywise absolute error: %g\n", err);

	// clean up
	algn_free(inv_expK_ref);
	algn_free(expK_ref);
	DeleteKineticExponential(&kinetic);
	DeleteSimulationParameters(&params);

	return (err < 1e-14 ? 0 : 1);
}
