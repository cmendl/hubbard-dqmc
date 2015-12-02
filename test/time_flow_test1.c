#include "time_flow.h"
#include "kinetic.h"
#include "util.h"
#include <mkl.h>
#include <math.h>
#include <stdio.h>


int TimeFlowTest1()
{
	// three-orbital simulation parameters
	sim_params_t params = { 0 };
	AllocateSimulationParameters(3, &params);

	// lattice dimensions
	params.Nx = 4;
	params.Ny = 5;

	int status;

	// read hopping parameters from disk
	status = ReadData("../test/time_flow_test1_taa.dat", params.t.aa, sizeof(double), params.Norb*params.Norb);	if (status < 0) { return status; }
	status = ReadData("../test/time_flow_test1_tab.dat", params.t.ab, sizeof(double), params.Norb*params.Norb);	if (status < 0) { return status; }
	status = ReadData("../test/time_flow_test1_tac.dat", params.t.ac, sizeof(double), params.Norb*params.Norb);	if (status < 0) { return status; }
	status = ReadData("../test/time_flow_test1_tad.dat", params.t.ad, sizeof(double), params.Norb*params.Norb);	if (status < 0) { return status; }
	status = ReadData("../test/time_flow_test1_tbc.dat", params.t.bc, sizeof(double), params.Norb*params.Norb);	if (status < 0) { return status; }

	// chemical potential
	params.mu = 5.0/7;

	// site energies
	params.eps[0] = -1.0/3;
	params.eps[1] =  2.0/11;
	params.eps[2] =  4.0/17;

	// imaginary-time step size
	params.dt = 1.0/8;

	// number of time steps
	params.L = 16;

	// largest number of B matrices multiplied together before performing a QR decomposition; must divide L
	params.prodBlen = 8;

	// lambda parameter (depends on U in Hamiltonian)
	const double lambda[] = { 0.75, 0.4, 1.0/9 };

	// total number of orbitals
	const int N = params.Nx*params.Ny * params.Norb;

	// calculate matrix exponential of the kinetic nearest neighbor hopping matrix
	kinetic_t kinetic;
	RectangularKineticExponential(&params, &kinetic);

	// load Hubbard-Stratonovich field from disk
	spin_field_t *s = (spin_field_t *)MKL_malloc(N*params.L *sizeof(spin_field_t), MEM_DATA_ALIGN);
	status = ReadData("../test/time_flow_test1_HS.dat", s, sizeof(spin_field_t), N*params.L); if (status != 0) { return status; }

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

	// allocate structures for output, represented as Q * diag(d) * T
	double *Q   = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	double *tau = (double *)MKL_malloc(N   * sizeof(double), MEM_DATA_ALIGN);	// scalar factors of the elementary reflectors for the matrix Q
	double *d   = (double *)MKL_malloc(N   * sizeof(double), MEM_DATA_ALIGN);
	double *T   = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);

	// calculate the imaginary-time flow map
	printf("Calculating imaginary-time flow map on a %i x %i lattice with %i orbitals per unit cell...\n", params.Nx, params.Ny, params.Norb);
	TimeFlowMap(&tsm, 0, Q, tau, d, T);

	// form A = Q * diag(d) * T (only for comparison with reference)
	for (i = 0; i < N; i++)
	{
		cblas_dscal(N, d[i], &T[i], N);
	}
	// multiply by Q from the left
	LAPACKE_dormqr(LAPACK_COL_MAJOR, 'L', 'N', N, N, N, Q, N, tau, T, N);

	// load reference data from disk
	double *A_ref = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	ReadData("../test/time_flow_test1_A.dat", A_ref, sizeof(double), N*N);

	// entrywise relative error
	double err_rel = 0;
	for (i = 0; i < N*N; i++)
	{
		err_rel = fmax(err_rel, fabs((T[i] - A_ref[i])/A_ref[i]));
	}
	printf("Largest entrywise relative error: %g\n", err_rel);

	// entrywise absolute error
	double err_abs = UniformDistance(N*N, T, A_ref);
	printf("Largest entrywise absolute error: %g\n", err_abs);

	// clean up
	MKL_free(A_ref);
	MKL_free(T);
	MKL_free(d);
	MKL_free(tau);
	MKL_free(Q);
	DeleteTimeStepMatrices(&tsm);
	MKL_free(expV[1]);
	MKL_free(expV[0]);
	MKL_free(s);
	DeleteKineticExponential(&kinetic);
	DeleteSimulationParameters(&params);

	return (err_rel < 4e-11 && err_abs < 1e-8 ? 0 : 1);
}
