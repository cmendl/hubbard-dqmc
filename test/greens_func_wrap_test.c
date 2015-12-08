#include "greens_func.h"
#include "kinetic.h"
#include "util.h"
#include <mkl.h>
#include <math.h>
#include <stdio.h>


extern void ComputeTimeStepMatrix(const kinetic_t *restrict kinetic, const double *const expV[2], const spin_field_t *restrict s, double *restrict B);

extern void ComputeInverseTimeStepMatrix(const kinetic_t *restrict kinetic, const double *const expV[2], const spin_field_t *restrict s, double *restrict invB);


int GreensFuncWrapTest()
{
	sim_params_t params = { 0 };
	AllocateSimulationParameters(1, &params);

	// lattice field dimensions
	params.Nx = 4;
	params.Ny = 6;

	// total number of lattice sites
	const int N = params.Norb * params.Nx*params.Ny;

	// imaginary-time step size
	params.dt = 1.0/8;

	// hopping parameters
	params.t.aa[0] = 0.0;
	params.t.ab[0] = 1.0;
	params.t.ac[0] = 1.0;
	params.t.ad[0] = 0.0;
	params.t.bc[0] = 0.0;

	// chemical potential
	params.mu = -1.0/5.0;
	params.eps[0] = 0;

	// lambda parameter (depends on U in Hamiltonian)
	const double lambda = 0.75;
	const double expV0 = exp(-lambda);
	const double expV1 = exp( lambda);
	const double *expV[2] = {&expV0, &expV1};

	// calculate matrix exponential of the kinetic nearest neighbor hopping matrix
	kinetic_t kinetic;
	RectangularKineticExponential(&params, &kinetic);

	// Hubbard-Stratonovich field
	// TODO: replace with malloc and load this from a file
	const spin_field_t s[4 * 6] = { 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1 };

	// calculate B and B^{-1}
	double *B    = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	double *invB = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	ComputeTimeStepMatrix(&kinetic, expV, s, B);
	ComputeInverseTimeStepMatrix(&kinetic, expV, s, invB);

	// load initial "Green's function" matrix from disk
	double *G = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	ReadData("../test/greens_func_wrap_test_G0.dat", G, sizeof(double), N*N);

	// perform a time slice wrap of the Green's function matrix
	printf("Performing a time slice wrap of the Green's function matrix on a %i x %i lattice...\n", params.Nx, params.Ny);
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
	DeleteSimulationParameters(&params);

	return (err_rel < 4e-12 && err_abs < 2e-14 ? 0 : 1);
}
