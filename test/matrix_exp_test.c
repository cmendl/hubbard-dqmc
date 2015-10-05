#include "linalg.h"
#include <math.h>
#include <stdio.h>


int MatrixExpTest()
{
	#define n 5
	const double A[n*n] = {
		 2,   -0.5,   1,    2.5, -1.5,
		-0.5,  3,     0.5,  2.5,  1.5,
		 1,    0.5,  -3,   -1,   -1,
		 2.5,  2.5,  -1,   -2.5,  2,
		-1.5,  1.5,  -1,    2,    1 };

	// reference entries of exp(A)
	const double expA_ref[n*n] = {
		 25.955160669119373,  -9.770400042675018, 4.871553760371625,   2.7134175545724095, -15.334376171601765,
		 -9.770400042675028, 105.8959791472407,  -8.733225951221918,  50.05805070380273,    66.66275673512335,
		  4.871553760371625,  -8.733225951221916, 1.6299718135475225, -3.2373011106885214,  -7.6198789670106155,
		  2.713417554572406,  50.05805070380271, -3.2373011106885223, 26.201262114829625,   29.683741571682877,
		-15.334376171601768,  66.66275673512334, -7.6198789670106155, 29.683741571682877,   47.30452009463436 };

	printf("Calculating matrix exponential of a %i x %i matrix...\n", n, n);
	double expA[n*n];
	MatrixExp(n, A, expA);

	// compare with reference
	double err = 0;
	int i;
	for (i = 0; i < n*n; i++)
	{
		err = fmax(err, fabs((expA[i] - expA_ref[i])/expA_ref[i]));
	}
	printf("Largest entrywise relative error: %g\n", err);

	return (err < 4e-15 ? 0 : 1);
}
