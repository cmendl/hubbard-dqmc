#include "greens_func.h"
#include "linalg.h"
#include "util.h"
#include <math.h>
#include <stdio.h>


int GreenUnequalTimeTest()
{
	int status;

	const int N = 5;            // block size
	const int L = 6;            // number of time steps
	const int prodBlen = 3;     // number of B matrices pre-multiplied together

	printf("Computing the unequal time Green's function for N = %i and L = %i...\n", N, L);

	// allocate time step matrices
	time_step_matrices_t tsm;
	AllocateTimeStepMatrices(N, L, prodBlen, &tsm);

	// load B matrices and their inverses from disk
	int l;
	for (l = 0; l < L; l++)
	{
		char path[1024];

		sprintf(path, "../test/green_unequal_time_test_B%i.dat", l);
		status = ReadData(path, tsm.B[l], sizeof(double), N*N);
		if (status != 0) { return status; }

		sprintf(path, "../test/green_unequal_time_test_invB%i.dat", l);
		status = ReadData(path, tsm.invB[l], sizeof(double), N*N);
		if (status != 0) { return status; }
	}

	// form products of the B matrices
	for (l = 0; l < tsm.numBprod; l++)
	{
		MatrixProductSequence(N, prodBlen, (tsm.B + l*prodBlen), tsm.Bprod[l]);
	}

	// temporary H matrix
	double *H = (double *)algn_malloc(N*N*tsm.numBprod*tsm.numBprod * sizeof(double));

	// compute unequal time Green's functions
	double *Gtau0 = (double *)algn_malloc(L*N*N * sizeof(double));
	double *G0tau = (double *)algn_malloc(L*N*N * sizeof(double));
	double *Geqlt = (double *)algn_malloc(L*N*N * sizeof(double));
	ComputeUnequalTimeGreensFunction(N, L, &tsm, H, Gtau0, G0tau, Geqlt);

	// load reference data from disk
	double *Gtau0_ref = (double *)algn_malloc(L*N*N * sizeof(double));
	double *G0tau_ref = (double *)algn_malloc(L*N*N * sizeof(double));
	double *Geqlt_ref = (double *)algn_malloc(L*N*N * sizeof(double));
	status = ReadData("../test/green_unequal_time_test_Gtau0.dat", Gtau0_ref, sizeof(double), L*N*N); if (status != 0) { return status; }
	status = ReadData("../test/green_unequal_time_test_G0tau.dat", G0tau_ref, sizeof(double), L*N*N); if (status != 0) { return status; }
	status = ReadData("../test/green_unequal_time_test_Geqlt.dat", Geqlt_ref, sizeof(double), L*N*N); if (status != 0) { return status; }

	// entrywise absolute error of matrix entries
	double err = fmax(fmax(
		UniformDistance(L*N*N, Gtau0, Gtau0_ref),
		UniformDistance(L*N*N, G0tau, G0tau_ref)),
		UniformDistance(L*N*N, Geqlt, Geqlt_ref));
	printf("Largest entrywise absolute error: %g\n", err);

	// clean up
	algn_free(Geqlt_ref);
	algn_free(G0tau_ref);
	algn_free(Gtau0_ref);
	algn_free(Geqlt);
	algn_free(G0tau);
	algn_free(Gtau0);
	algn_free(H);
	DeleteTimeStepMatrices(&tsm);

	return (err < 4e-14 ? 0 : 1);
}
