#include "kinetic.h"
#include "util.h"
#include <mkl.h>
#include <math.h>
#include <stdio.h>


int KineticTest()
{
	// imaginary-time step size
	const double dt = 1.0/7;

	#define Nx 4
	#define Ny 6
	// total number of lattice sites
	#define N  (Nx * Ny)

	// calculate matrix exponential of the kinetic nearest neighbor hopping matrix
	printf("Calculating matrix exponential of the kinetic nearest neighbor hopping matrix on a %i x %i lattice...\n", Nx, Ny);
	kinetic_t kinetic;
	NearestNeighborKineticExponential(Nx, Ny, 0.0, dt, &kinetic);

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

	return (err < 2e-15 ? 0 : 1);
}
