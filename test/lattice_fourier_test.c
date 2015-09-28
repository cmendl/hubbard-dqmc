#include "lattice_fourier.h"
#include "util.h"
#include <mkl.h>
#include <math.h>
#include <stdio.h>


//________________________________________________________________________________________________________________________
///
/// \brief Complex conjugate
///
static inline MKL_Complex16 ComplexConjugate(const MKL_Complex16 z)
{
	MKL_Complex16 ret;
	ret.real =  z.real;
	ret.imag = -z.imag;

	return ret;
}


//________________________________________________________________________________________________________________________
///
/// \brief Subtract two complex numbers
///
static inline MKL_Complex16 ComplexSubtract(const MKL_Complex16 z, const MKL_Complex16 w)
{
	MKL_Complex16 ret;
	ret.real = z.real - w.real;
	ret.imag = z.imag - w.imag;

	return ret;
}


//________________________________________________________________________________________________________________________
///
/// \brief Absolute value |z| of a complex number z
///
static inline double ComplexAbs(const MKL_Complex16 z)
{
	if (fabs(z.real) < fabs(z.imag))
	{
		return fabs(z.imag)*sqrt(1.0 + square(z.real/z.imag));
	}
	else	// fabs(z.imag) <= fabs(z.real)
	{
		if (z.real != 0.0)
		{
			return fabs(z.real)*sqrt(1.0 + square(z.imag/z.real));
		}
		else
		{
			return 0.0;
		}
	}
}


int LatticeFourierTestDim(const int Nx, const int Ny)
{
	// total number of lattice sites
	const int N = Nx * Ny;

	lattice_fourier_desc_t desc;
	CreateLatticeFourierDescriptor(Nx, Ny, &desc);

	// load input vector from disk
	double *v = (double *)MKL_malloc(N * sizeof(double), MEM_DATA_ALIGN);
	char fn[1024];
	sprintf(fn, "../test/lattice_fourier_test_%ix%i_v.dat", Nx, Ny);
	ReadData(fn, v, sizeof(double), N);

	// compute forward Fourier transform
	printf("Computing the forward Fourier transform of a real vector on a %i x %i lattice...\n", Nx, Ny);
	MKL_Complex16 *vF = (MKL_Complex16 *)MKL_malloc(Ny*(Nx/2+1) * sizeof(MKL_Complex16), MEM_DATA_ALIGN);
	FourierTransformLatticeForward(&desc, v, vF);

	// Fourier transform back
	printf("Computing the backward Fourier transform of a complex vector on a %i x %i lattice...\n", Nx, Ny);
	double *v1 = (double *)MKL_malloc(N * sizeof(double), MEM_DATA_ALIGN);
	FourierTransformLatticeBackward(&desc, vF, v1);
	// normalization
	int i;
	for (i = 0; i < N; i++)
	{
		v1[i] *= desc.nfac;
	}

	// load reference data from disk
	MKL_Complex16 *vF_ref = (MKL_Complex16 *)MKL_malloc(N * sizeof(MKL_Complex16), MEM_DATA_ALIGN);
	sprintf(fn, "../test/lattice_fourier_test_%ix%i_vF.dat", Nx, Ny);
	ReadData(fn, vF_ref, sizeof(MKL_Complex16), N);

	// compare 'vF' with reference
	double errF = 0;
	int j;
	for (j = 0; j < Ny; j++)
	{
		for (i = 0; i <= Nx/2; i++)
		{
			errF = fmax(errF, ComplexAbs(ComplexSubtract(vF[i + j*(Nx/2+1)], vF_ref[i + j*Nx])));
		}

		// effectively flip the sign of j
		const int j_neg = (Ny - j) % Ny;

		for (i = Nx/2 + 1; i < Nx; i++)
		{
			errF = fmax(errF, ComplexAbs(ComplexSubtract(ComplexConjugate(vF[(Nx-i) + j_neg*(Nx/2+1)]), vF_ref[i + j*Nx])));
		}
	}
	printf("Largest entrywise absolute error of forward Fourier transform: %g\n", errF);

	// compare 'v1' with 'v'
	double err1 = UniformDistance(N, v1, v);
	printf("Largest entrywise absolute error after forward and backward Fourier transform: %g\n", err1);

	// clean up
	MKL_free(vF_ref);
	MKL_free(v1);
	MKL_free(vF);
	MKL_free(v);
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
