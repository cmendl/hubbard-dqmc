#ifndef LATTICE_FOURIER_H
#define LATTICE_FOURIER_H

#ifdef USE_MKL

#include <complex.h>
#include <mkl_dfti.h>


//________________________________________________________________________________________________________________________
///
/// \brief FFT descriptors for performing forward and backward Fourier transforms on the rectangular lattice
///
typedef struct
{
	DFTI_DESCRIPTOR_HANDLE hand_forw;   //!< forward  Fourier transform descriptor
	DFTI_DESCRIPTOR_HANDLE hand_back;   //!< backward Fourier transform descriptor

	double nfac;                        //!< normalization factor

	int Nx;                             //!< x-dimension of the square lattice
	int Ny;                             //!< y-dimension of the square lattice
}
lattice_fourier_desc_t;


MKL_LONG CreateLatticeFourierDescriptor(const int Nx, const int Ny, lattice_fourier_desc_t *desc);

void DeleteLatticeFourierDescriptor(lattice_fourier_desc_t *desc);


//________________________________________________________________________________________________________________________
///
/// \brief Perform a forward Fourier transformation on the lattice
///
inline void FourierTransformLatticeForward(lattice_fourier_desc_t *desc, const double *v, double complex *vF)
{
	DftiComputeForward(desc->hand_forw, (double *)v, vF);
}


//________________________________________________________________________________________________________________________
///
/// \brief Perform a reverse Fourier transformation on the lattice
///
inline void FourierTransformLatticeBackward(lattice_fourier_desc_t *desc, const double complex *vF, double *v)
{
	DftiComputeBackward(desc->hand_back, (double *)vF, v);
}


#endif

#endif
