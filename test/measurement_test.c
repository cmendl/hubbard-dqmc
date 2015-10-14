#include "measurement.h"
#include "util.h"
#include <mkl.h>
#include <math.h>
#include <stdio.h>


int MeasurementTest()
{
	int status;

	const int N = 5;	// block size
	const int L = 4;	// number of time steps

	printf("Performing an unequal time measurement for N = %i and L = %i...\n", N, L);

	measurement_data_unequal_time_t meas_data;
	AllocateUnequalTimeMeasurementData(N, L, &meas_data);

	// load B matrices from disk
	double **Bu = (double **)MKL_malloc(L * sizeof(double *), MEM_DATA_ALIGN);
	double **Bd = (double **)MKL_malloc(L * sizeof(double *), MEM_DATA_ALIGN);
	int l;
	for (l = 0; l < L; l++)
	{
		Bu[l] = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
		Bd[l] = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);

		char path[1024];
		sprintf(path, "../test/measurement_test_Bu%i.dat", l);  status = ReadData(path, Bu[l], sizeof(double), N*N); if (status != 0) { return status; }
		sprintf(path, "../test/measurement_test_Bd%i.dat", l);  status = ReadData(path, Bd[l], sizeof(double), N*N); if (status != 0) { return status; }
	}

	// accumulate an unequal time measurement (compute Green's functions)
	AccumulateUnequalTimeMeasurement(1, Bu, Bd, &meas_data);

	// load reference data from disk
	double *Gu_ref = (double *)MKL_malloc(L*N*N * sizeof(double), MEM_DATA_ALIGN);
	double *Gd_ref = (double *)MKL_malloc(L*N*N * sizeof(double), MEM_DATA_ALIGN);
	status = ReadData("../test/measurement_test_Gu.dat", Gu_ref, sizeof(double), L*N*N); if (status != 0) { return status; }
	status = ReadData("../test/measurement_test_Gd.dat", Gd_ref, sizeof(double), L*N*N); if (status != 0) { return status; }

	// entrywise absolute error of matrix entries
	double err = fmax(
		UniformDistance(L*N*N, meas_data.Gu_ut, Gu_ref),
		UniformDistance(L*N*N, meas_data.Gd_ut, Gd_ref));
	printf("Largest entrywise absolute error: %g\n", err);

	// clean up
	MKL_free(Gd_ref);
	MKL_free(Gu_ref);
	for (l = 0; l < L; l++)
	{
		MKL_free(Bd[l]);
		MKL_free(Bu[l]);
	}
	MKL_free(Bd);
	MKL_free(Bu);
	DeleteUnequalTimeMeasurementData(&meas_data);

	return (err < 1e-14 ? 0 : 1);
}
