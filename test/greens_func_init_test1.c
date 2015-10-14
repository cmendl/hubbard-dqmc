#include "greens_func.h"
#include "kinetic.h"
#include "util.h"
#include <mkl.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>


int GreensFuncInitTest1()
{
	int i;

	// lattice field dimensions
	#define Nx 4
	#define Ny 6
	// total number of lattice sites
	#define N  (Nx * Ny)

	// imaginary-time step size
	const double dt = 1.0/8;

	// number of time steps
	#define L 16

	// largest number of B_l matrices multiplied together before performing a QR decomposition; must divide L
	const int prodBlen = 4;

	const double lambda = 0.75;
	const double expV[2] = {
		exp(-lambda),
		exp( lambda)
	};

	// calculate matrix exponential of the kinetic nearest neighbor hopping matrix
	kinetic_t kinetic;
	NearestNeighborKineticExponential(Nx, Ny, 0.0, dt, &kinetic);

	// Hubbard-Stratonovich field
	const spin_field_t s[L*N] = {
		1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1,
		0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0,
		0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1,
		1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1,
		0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0,
		0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0,
		1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1,
		0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0,
		1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1,
		1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0,
		1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0,
		1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0,
		0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1,
		1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1,
		1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0,
		1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1
	};

	// allocate and initialize time step matrices
	time_step_matrices_t tsm;
	AllocateTimeStepMatrices(N, L, prodBlen, &tsm);
	InitTimeStepMatrices(&kinetic, expV, s, &tsm);

	// construct the Green's function matrix
	printf("Constructing Green's function for a %i x %i lattice at beta = %g...\n", Nx, Ny, L*dt);
	greens_func_t G;
	AllocateGreensFunction(N, &G);
	GreenConstruct(&tsm, 0, &G);

	// load reference data from disk
	double Gmat_ref[N*N];
	double detG_ref;
	ReadData("../test/greens_func_init_test1_G.dat", Gmat_ref, sizeof(double), N*N);
	ReadData("../test/greens_func_init_test1_detG.dat", &detG_ref, sizeof(double), 1);

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
	DeleteGreensFunction(&G);
	DeleteTimeStepMatrices(&tsm);
	DeleteKineticExponential(&kinetic);

	return (err_rel < 4e-12 && err_abs < 2e-14 && err_det < 5e-14 ? 0 : 1);
}
