#include "time_flow.h"
#include "kinetic.h"
#include "util.h"
#include <mkl.h>
#include <math.h>
#include <stdio.h>


int TimeFlowTest3()
{
	// two-orbital simulation parameters
	sim_params_t params = { 0 };
	AllocateSimulationParameters(2, &params);

	// lattice field dimensions
	params.Nx = 4;
	params.Ny = 6;

	// total number of orbitals
	const int N = params.Norb * params.Nx*params.Ny;

	// imaginary-time step size
	params.dt = 1.0/8;

	// t (nearest neighbor) hopping parameter
	params.t.aa[1] = -2.0/3;
	// a -> b
	params.t.ab[0] =  4.0/9;
	params.t.ab[1] =  3.0/7;
	params.t.ab[2] =  3.0/16;
	params.t.ab[3] = -2.0/19;
	// a -> c
	params.t.ac[0] =  1.0/5;
	params.t.ac[1] =  2.0/9;
	params.t.ac[2] = -1.0/9;
	params.t.ac[3] = -3.0/7;
	// a -> d
	params.t.ad[0] =  1.0/13;
	params.t.ad[1] = -1.0/10;
	params.t.ad[2] =  4.0/5;
	params.t.ad[3] =  2.0/17;
	// b -> c
	params.t.bc[0] =  3.0/19;
	params.t.bc[1] =  1.0/6;
	params.t.bc[2] = -3.0/11;
	params.t.bc[3] =  1.0/7;

	// chemical potential
	params.mu = 6.0/7;
	// site energies
	params.eps[0] =  1.0/5;
	params.eps[1] = -2.0/11;

	// electron-phonon interaction strength
	params.phonon_params.g[0] = 1.0/7;
	params.phonon_params.g[1] = 3.0/13;

	// lambda parameter (depends on U in Hamiltonian)
	const double lambda[] = { 0.75, 5.0/23 };

	// number of time steps
	params.L = 16;

	// largest number of B matrices multiplied together before performing a QR decomposition; must divide L
	params.prodBlen = 8;

	// calculate matrix exponential of the kinetic nearest neighbor hopping matrix
	kinetic_t kinetic;
	RectangularKineticExponential(&params, &kinetic);

	int status;

	// load Hubbard-Stratonovich field from disk
	spin_field_t *s = (spin_field_t *)MKL_malloc(N*params.L *sizeof(spin_field_t), MEM_DATA_ALIGN);
	status = ReadData("../test/time_flow_test3_HS.dat", s, sizeof(spin_field_t), N*params.L);
	if (status != 0) { return status; }

	// interaction potential
	double *expV[2];
	expV[0] = (double *)MKL_malloc(params.Norb * sizeof(double), MEM_DATA_ALIGN);
	expV[1] = (double *)MKL_malloc(params.Norb * sizeof(double), MEM_DATA_ALIGN);
	int i;
	for (i = 0; i < params.Norb; i++)
	{
		expV[0][i] = exp(-lambda[i]);
		expV[1][i] = exp( lambda[i]);
	}

	// phonon field
	double *X    = (double *)MKL_malloc(params.L*N * sizeof(double), MEM_DATA_ALIGN);
	double *expX = (double *)MKL_malloc(params.L*N * sizeof(double), MEM_DATA_ALIGN);
	status = ReadData("../test/time_flow_test3_X.dat", X, sizeof(double), params.L*N);
	if (status != 0) { return status; }
	const int Ncell = params.Nx * params.Ny;
	int l;
	for (l = 0; l < params.L; l++)
	{
		for (i = 0; i < N; i++)
		{
			const int o = i / Ncell;	// orbital number
			expX[i + l*N] = exp(-params.dt*params.phonon_params.g[o] * X[i + l*N]);
		}
	}

	// allocate and initialize time step matrices
	time_step_matrices_t tsm;
	AllocateTimeStepMatrices(N, params.L, params.prodBlen, &tsm);
	InitPhononTimeStepMatrices(&kinetic, expV, s, expX, &tsm);

	// allocate structures for output, represented as Q * diag(d) * T
	double *Q   = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	double *tau = (double *)MKL_malloc(N   * sizeof(double), MEM_DATA_ALIGN);	// scalar factors of the elementary reflectors for the matrix Q
	double *d   = (double *)MKL_malloc(N   * sizeof(double), MEM_DATA_ALIGN);
	double *T   = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);

	// calculate the imaginary-time flow map
	printf("Calculating imaginary-time flow map including phonons on a %i x %i lattice...\n", params.Nx, params.Ny);
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
	ReadData("../test/time_flow_test3_A.dat", A_ref, sizeof(double), N*N);

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
	MKL_free(expX);
	MKL_free(X);
	DeleteKineticExponential(&kinetic);

	return (err_rel < 1e-10 && err_abs < 4e-10 ? 0 : 1);
}
