#include "greens_func.h"
#include "kinetic.h"
#include "util.h"
#include <mkl.h>
#include <math.h>
#include <stdio.h>


void ComputeTimeStepMatrix(const kinetic_t *restrict kinetic, const double expV[2], const spin_field_t *restrict s, double *restrict B);

void ComputeInverseTimeStepMatrix(const kinetic_t *restrict kinetic, const double expV[2], const spin_field_t *restrict s, double *restrict invB);


int GreensFuncWrapTest()
{
	// lattice field dimensions
	#define Nx 4
	#define Ny 6
	// total number of lattice sites
	#define N  (Nx * Ny)

	// imaginary-time step size
	const double dt = 1.0/8;

	// chemical potential
	const double mu = -1.0/5;

	// lambda parameter (depends on U in Hamiltonian)
	const double lambda = 0.75;
	const double expV[2] = {
		exp(-lambda),
		exp( lambda)
	};

	// calculate matrix exponential of the kinetic nearest neighbor hopping matrix
	kinetic_t kinetic;
	NearestNeighborKineticExponential(Nx, Ny, mu, dt, &kinetic);

	// Hubbard-Stratonovich field
	const spin_field_t s[N] = { 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1 };

	// calculate B and B^{-1}
	double *B    = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	double *invB = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	ComputeTimeStepMatrix(&kinetic, expV, s, B);
	ComputeInverseTimeStepMatrix(&kinetic, expV, s, invB);

	// load initial "Green's function" matrix from disk
	double *G = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	ReadData("../test/greens_func_wrap_test_G0.dat", G, sizeof(double), N*N);

	// perform a time slice wrap of the Green's function matrix
	printf("Performing a time slice wrap of the Green's function matrix on a %i x %i lattice...\n", Nx, Ny);
	GreenTimeSliceWrap(N, B, invB, G);

	// load reference data from disk
	double G_ref[N*N];
	ReadData("../test/greens_func_wrap_test_G1.dat", G_ref, sizeof(double), N*N);

	// entrywise relative error
	double err_rel = 0;
	int i;
	for (i = 0; i < N*N; i++)
	{
		err_rel = fmax(err_rel, fabs((G[i] - G_ref[i])/G_ref[i]));
	}
	printf("Largest entrywise relative error: %g\n", err_rel);

	// entrywise absolute error
	double err_abs = UniformDistance(N*N, G, G_ref);
	printf("Largest entrywise absolute error: %g\n", err_abs);

	// clean up
	MKL_free(G);
	MKL_free(invB);
	MKL_free(B);
	DeleteKineticExponential(&kinetic);

	return (err_rel < 4e-12 && err_abs < 2e-14 ? 0 : 1);
}
