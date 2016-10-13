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
	// three-orbital simulation parameters
	sim_params_t params = { 0 };
	AllocateSimulationParameters(3, &params);

	// lattice field dimensions
	params.Nx = 4;
	params.Ny = 6;

	// total number of lattice sites
	const int N = params.Norb * params.Nx*params.Ny;

	// imaginary-time step size
	params.dt = 1.0/8.0;

	int status;

	// read hopping parameters from disk
	status = ReadData("../test/greens_func_init_test1_taa.dat", params.t.aa, sizeof(double), params.Norb*params.Norb);	if (status < 0) { return status; }
	status = ReadData("../test/greens_func_init_test1_tab.dat", params.t.ab, sizeof(double), params.Norb*params.Norb);	if (status < 0) { return status; }
	status = ReadData("../test/greens_func_init_test1_tac.dat", params.t.ac, sizeof(double), params.Norb*params.Norb);	if (status < 0) { return status; }
	status = ReadData("../test/greens_func_init_test1_tad.dat", params.t.ad, sizeof(double), params.Norb*params.Norb);	if (status < 0) { return status; }
	status = ReadData("../test/greens_func_init_test1_tbc.dat", params.t.bc, sizeof(double), params.Norb*params.Norb);	if (status < 0) { return status; }

	// chemical potential
	params.mu = 5.0/6;

	// site energies
	params.eps[0] = -1.0/13;
	params.eps[1] =  5.0/17;
	params.eps[2] =  4.0/11;

	// number of time steps
	params.L = 16;

	// largest number of B_l matrices multiplied together before performing a QR decomposition; must divide L
	params.prodBlen = 4;

	// lambda parameter (depends on U in Hamiltonian)
	const double lambda[] = { 0.5, 2.0/9, 0.2 };

	// calculate matrix exponential of the kinetic hopping matrix
	kinetic_t kinetic;
	RectangularKineticExponential(&params, &kinetic);

	// load Hubbard-Stratonovich field from disk
	spin_field_t *s = (spin_field_t *)MKL_malloc(N*params.L *sizeof(spin_field_t), MEM_DATA_ALIGN);
	status = ReadData("../test/greens_func_init_test1_HS.dat", s, sizeof(spin_field_t), N*params.L); if (status != 0) { return status; }

	double *expV[2];
	expV[0] = (double *)MKL_malloc(params.Norb * sizeof(double), MEM_DATA_ALIGN);
	expV[1] = (double *)MKL_malloc(params.Norb * sizeof(double), MEM_DATA_ALIGN);
	int i;
	for (i = 0; i < params.Norb; i++)
	{
		expV[0][i] = exp(-lambda[i]);
		expV[1][i] = exp( lambda[i]);
	}

	// allocate and initialize time step matrices
	time_step_matrices_t tsm;
	AllocateTimeStepMatrices(N, params.L, params.prodBlen, &tsm);
	InitTimeStepMatrices(&kinetic, expV, s, &tsm);

	// construct the Green's function matrix
	printf("Constructing Green's function for a %i x %i lattice with %i orbitals per unit cell at beta = %g...\n", params.Nx, params.Ny, params.Norb, params.L*params.dt);
	greens_func_t G;
	AllocateGreensFunction(N, &G);
	GreenConstruct(&tsm, 0, &G);

	// load reference data from disk
	double *Gmat_ref = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
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
	MKL_free(Gmat_ref);
	DeleteGreensFunction(&G);
	DeleteTimeStepMatrices(&tsm);
	DeleteKineticExponential(&kinetic);
	MKL_free(expV[1]);
	MKL_free(expV[0]);
	MKL_free(s);
	DeleteSimulationParameters(&params);

	return (err_rel < 2e-10 && err_abs < 1e-14 && err_det < 2e-13 ? 0 : 1);
}
