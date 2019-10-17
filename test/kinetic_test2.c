#include "kinetic.h"
#include "util.h"
#include <math.h>
#include <stdio.h>


int KineticTest2()
{
	// three-orbital simulation parameters
	sim_params_t params = { 0 };
	AllocateSimulationParameters(3, &params);

	// lattice dimension
	params.Nx = 5;
	params.Ny = 4;

	int status;

	// read hopping parameters from disk
	status = ReadData("../test/kinetic_test2_taa.dat", params.t.aa, sizeof(double), params.Norb*params.Norb);	if (status < 0) { return status; }
	status = ReadData("../test/kinetic_test2_tab.dat", params.t.ab, sizeof(double), params.Norb*params.Norb);	if (status < 0) { return status; }
	status = ReadData("../test/kinetic_test2_tac.dat", params.t.ac, sizeof(double), params.Norb*params.Norb);	if (status < 0) { return status; }
	status = ReadData("../test/kinetic_test2_tad.dat", params.t.ad, sizeof(double), params.Norb*params.Norb);	if (status < 0) { return status; }
	status = ReadData("../test/kinetic_test2_tbc.dat", params.t.bc, sizeof(double), params.Norb*params.Norb);	if (status < 0) { return status; }

	// chemical potential
	params.mu = 5.0/7;

	// site energies
	params.eps[0] = -2.0/11;
	params.eps[1] =  1.0/3;
	params.eps[2] =  4.0/13;

	// imaginary-time step size
	params.dt = 1.0/9;

	// set some entries to pass parameter validation (actual values not used here)
	params.L        = 1;
	params.prodBlen = 1;
	params.nwraps   = 1;
	params.itime    = 1;

	status = ValidateSimulationParameters(&params);
	if (status < 0) { return status; }

	// calculate matrix exponential of the kinetic hopping matrix
	printf("Calculating matrix exponential of the kinetic hopping matrix on a %i x %i lattice with %i orbitals...\n", params.Nx, params.Ny, params.Norb);
	kinetic_t kinetic;
	RectangularKineticExponential(&params, &kinetic);

	// load reference data from disk
	const int N = kinetic.Ncell * kinetic.Norb;
	double *expK_ref     = (double *)algn_malloc(N*N * sizeof(double));
	double *inv_expK_ref = (double *)algn_malloc(N*N * sizeof(double));
	status = ReadData("../test/kinetic_test2_expK.dat",    expK_ref,     sizeof(double), N*N);	if (status < 0) { return status; }
	status = ReadData("../test/kinetic_test2_invexpK.dat", inv_expK_ref, sizeof(double), N*N);	if (status < 0) { return status; }

	// compare with reference
	double err = fmax(
		UniformDistance(N*N, kinetic.expK,         expK_ref),
		UniformDistance(N*N, kinetic.inv_expK, inv_expK_ref));
	printf("Largest entrywise absolute error: %g\n", err);

	// clean up
	algn_free(inv_expK_ref);
	algn_free(expK_ref);
	DeleteKineticExponential(&kinetic);
	DeleteSimulationParameters(&params);

	return (err < 1e-14 ? 0 : 1);
}
