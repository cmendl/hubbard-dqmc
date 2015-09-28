#include "greens_func.h"
#include "util.h"
#include <math.h>
#include <stdio.h>


int GreensFuncFlipTest()
{
	// total number of lattice sites
	#define N 24

	// nonzero entry of the Delta matrix
	const double delta = 1.0/7;

	// load initial "Green's function" matrix from disk
	double G[N*N];
	ReadData("../test/greens_func_flip_test_G0.dat", G, sizeof(double), N*N);

	// update Green's function matrix after a spin flip
	printf("Updating Green's function matrix after a spin flip...\n");
	GreenShermanMorrisonUpdate(delta, N, 11, G);

	// load reference data from disk
	double G_ref[N*N];
	ReadData("../test/greens_func_flip_test_G1.dat", G_ref, sizeof(double), N*N);

	// entrywise absolute error
	double err = UniformDistance(N*N, G, G_ref);
	printf("Largest entrywise absolute error: %g\n", err);

	return (err < 2e-15 ? 0 : 1);
}
