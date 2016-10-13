#include "greens_func.h"
#include "kinetic.h"
#include "util.h"
#include <mkl.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>


int GreensFuncInitTest2()
{
	int i;
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
	params.t.ab[0] = 1;
	params.t.ac[0] = 1;
	params.t.ad[0] = 0.2;
	params.t.bc[0] = 0.2;

	// chemical potential
	params.mu = -3.0/17;
	params.eps[0] = 0;

	// number of time steps
	params.L = 32;

	// largest number of B_l matrices multiplied together before performing a QR decomposition; must divide L
	params.prodBlen = 4;

	const double lambda = 0.75;
	double *expV[2] = {
		(double *)MKL_malloc(sizeof(double), MEM_DATA_ALIGN),
		(double *)MKL_malloc(sizeof(double), MEM_DATA_ALIGN)
	};
	expV[0][0] = exp(-lambda);
	expV[1][0] = exp( lambda);

	// calculate matrix exponential of the kinetic nearest neighbor hopping matrix
	kinetic_t kinetic;
	RectangularKineticExponential(&params, &kinetic);

	// Hubbard-Stratonovich field
	// TODO: replace with malloc and load this from a file
	const spin_field_t s[32 * 4 * 6] = {
		1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1,
		0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1,
		0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0,
		1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0,
		0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0,
		0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0,
		1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0,
		0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0,
		1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0,
		1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0,
		1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1,
		1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0,
		0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1,
		1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1,
		1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1,
		1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0,
		1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0,
		1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1,
		0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0,
		1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1,
		1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1,
		1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0,
		1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1,
		1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1,
		0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1,
		0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0,
		0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0,
		0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0,
		1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1,
		1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0,
		1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0,
		1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1
	};

	// allocate and initialize time step matrices
	time_step_matrices_t tsm;
	AllocateTimeStepMatrices(N, params.L, params.prodBlen, &tsm);
	InitTimeStepMatrices(&kinetic, expV, s, &tsm);

	// construct the Green's function matrix
	printf("Constructing Green's function for a %i x %i lattice at beta = %g...\n", params.Nx, params.Ny, params.L*params.dt);
	greens_func_t G;
	AllocateGreensFunction(N, &G);
	GreenConstruct(&tsm, 0, &G);

	// load reference data from disk
	double *Gmat_ref = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	double detG_ref;
	ReadData("../test/greens_func_init_test2_G.dat", Gmat_ref, sizeof(double), N*N);
	ReadData("../test/greens_func_init_test2_detG.dat", &detG_ref, sizeof(double), 1);

	// entrywise relative error of matrix entries
	double err_rel = 0;
	for (i = 0; i < N*N; i++)
	{
		err_rel = fmax(err_rel, fabs((G.mat[i] - Gmat_ref[i])/Gmat_ref[i]));
	}
	printf("Largest entrywise relative error: %g\n", err_rel);

	// entrywise absolute error of matrix entries
	double err_abs = UniformDistance(N*N, G.mat, Gmat_ref);
	printf("Largest entrywise absolute error: %g\n", err_abs);

	// relative error of determinant
	double detG = G.sgndet * exp(G.logdet);
	double err_det = fabs((detG - detG_ref) / detG_ref);
	printf("Relative determinant error: %g\n", err_det);

	// clean up
	MKL_free(Gmat_ref);
	DeleteGreensFunction(&G);
	DeleteTimeStepMatrices(&tsm);
	DeleteKineticExponential(&kinetic);
	MKL_free(expV[1]);
	MKL_free(expV[0]);
	DeleteSimulationParameters(&params);

	return (err_rel < 2e-11 && err_abs < 2e-14 && err_det < 5e-13 ? 0 : 1);
}
