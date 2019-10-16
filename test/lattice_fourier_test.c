#ifdef USE_MKL

#include "lattice_fourier.h"
#include "util.h"
#include <math.h>
#include <complex.h>
#include <stdio.h>


static int LatticeFourierTestDim(const int Nx, const int Ny)
{
	// total number of lattice sites
	const int N = Nx * Ny;

	lattice_fourier_desc_t desc;
	CreateLatticeFourierDescriptor(Nx, Ny, &desc);

	// load input vector from disk
	double *v = (double *)algn_malloc(N * sizeof(double));
	char fn[1024];
	sprintf(fn, "../test/lattice_fourier_test_%ix%i_v.dat", Nx, Ny);
	ReadData(fn, v, sizeof(double), N);

	// compute forward Fourier transform
	printf("Computing the forward Fourier transform of a real vector on a %i x %i lattice...\n", Nx, Ny);
	double complex *vF = (double complex *)algn_malloc(Ny*(Nx/2+1) * sizeof(double complex));
	FourierTransformLatticeForward(&desc, v, vF);

	// Fourier transform back
	printf("Computing the backward Fourier transform of a complex vector on a %i x %i lattice...\n", Nx, Ny);
	double *v1 = (double *)algn_malloc(N * sizeof(double));
	FourierTransformLatticeBackward(&desc, vF, v1);
	// normalization
	int i;
	for (i = 0; i < N; i++)
	{
		v1[i] *= desc.nfac;
	}

	// load reference data from disk
	double complex *vF_ref = (double complex *)algn_malloc(N * sizeof(double complex));
	sprintf(fn, "../test/lattice_fourier_test_%ix%i_vF.dat", Nx, Ny);
	ReadData(fn, vF_ref, sizeof(double complex), N);

	// compare 'vF' with reference
	double errF = 0;
	int j;
	for (j = 0; j < Ny; j++)
	{
		for (i = 0; i <= Nx/2; i++)
		{
			errF = fmax(errF, cabs(vF[i + j*(Nx/2+1)] - vF_ref[i + j*Nx]));
		}

		// effectively flip the sign of j
		const int j_neg = (Ny - j) % Ny;

		for (i = Nx/2 + 1; i < Nx; i++)
		{
			errF = fmax(errF, cabs(conj(vF[(Nx-i) + j_neg*(Nx/2+1)]) - vF_ref[i + j*Nx]));
		}
	}
	printf("Largest entrywise absolute error of forward Fourier transform: %g\n", errF);

	// compare 'v1' with 'v'
	double err1 = UniformDistance(N, v1, v);
	printf("Largest entrywise absolute error after forward and backward Fourier transform: %g\n", err1);

	// clean up
	algn_free(vF_ref);
	algn_free(v1);
	algn_free(vF);
	algn_free(v);
	DeleteLatticeFourierDescriptor(&desc);

	return (errF < 1e-15 && err1 < 5e-16 ? 0 : 1);
}


int LatticeFourierTest()
{
	const int hr0 = LatticeFourierTestDim(4, 6);
	const int hr1 = LatticeFourierTestDim(5, 6);

	if (hr0 < 0 || hr1 < 0)
	{
		return -1;
	}
	else
	{
		return hr0 + hr1;
	}
}


#else // USE_MKL


// skip test (current implementation requires MKL)

int LatticeFourierTest()
{
	return 0;
}


#endif
