#include "linalg.h"
#include "util.h"
#include <mkl.h>
#include <math.h>
#include <stdio.h>


int MatrixDetTest()
{
	#define num 4

	// read reference determinant values from disk
	double det_ref[num];
	int status = ReadData("../test/matrix_det_test_ref.dat", det_ref, sizeof(double), num);

	double err = 0;
	int i;
	for (i = 0; i < num; i++)
	{
		// matrix dimension
		const int n = i + 4;

		double *A = (double *)MKL_malloc(n*n * sizeof(double), MEM_DATA_ALIGN);
		char filename[1024];
		sprintf(filename, "../test/matrix_det_test_A%i.dat", i + 1);
		status = ReadData(filename, A, sizeof(double), n*n);
		if (status != 0) {
			return status;
		}

		printf("Calculating the determinant and determinant sign of a %i x %i matrix...\n", n, n);
		double det = Determinant(n, A);
		int sign = DeterminantSign(n, A);

		// compare with reference
		err = fmax(err, fabs(det - det_ref[i]));
		err = fmax(err, (double)(abs(sign - Signum(det_ref[i]))));

		// clean up
		MKL_free(A);
	}

	printf("Largest absolute error: %g\n", err);

	return (err < 5e-16 ? 0 : 1);
}
