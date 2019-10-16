#include "lattice_fourier.h"

#ifdef USE_MKL


//________________________________________________________________________________________________________________________
///
/// \brief Generate FFT descriptors for performing forward and backward Fourier transforms on the rectangular lattice
///
MKL_LONG CreateLatticeFourierDescriptor(const int Nx, const int Ny, lattice_fourier_desc_t *desc)
{
	desc->Nx = Nx;
	desc->Ny = Ny;

	// flip Nx <-> Ny since MKL FFT "halves" the last dimension

	MKL_LONG status = 0;

	// create, configure and commit DFTI descriptor for forward transform
	{
		MKL_LONG sizes[2];
		sizes[0] = Ny;
		sizes[1] = Nx;
		status = DftiCreateDescriptor(&desc->hand_forw, DFTI_DOUBLE, DFTI_REAL, 2, sizes);
		if (status != 0) { return status; }

		// set conjugate-even storage; DFTI_COMPLEX_COMPLEX setting overwrites DFTI_PACKED_FORMAT
		status = DftiSetValue(desc->hand_forw, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
		if (status != 0) { return status; }

		// set configuration: out-of-place
		status = DftiSetValue(desc->hand_forw, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
		if (status != 0) { return status; }

		// set real input strides
		MKL_LONG strides[3];
		strides[0] = 0;
		strides[1] = Nx;
		strides[2] = 1;
		status = DftiSetValue(desc->hand_forw, DFTI_INPUT_STRIDES, strides);
		if (status != 0) { return status; }

		// set complex output strides
		strides[0] = 0;
		strides[1] = Nx/2+1;
		strides[2] = 1;
		status = DftiSetValue(desc->hand_forw, DFTI_OUTPUT_STRIDES, strides);
		if (status != 0) { return status; }

		// commit the descriptor
		status = DftiCommitDescriptor(desc->hand_forw);
		if (status != 0) { return status; }
	}

	// create, configure and commit DFTI descriptor for backward transform
	{
		MKL_LONG sizes[2];
		sizes[0] = Ny;
		sizes[1] = Nx;
		status = DftiCreateDescriptor(&desc->hand_back, DFTI_DOUBLE, DFTI_REAL, 2, sizes);
		if (status != 0) { return status; }

		// set conjugate-even storage; DFTI_COMPLEX_COMPLEX setting overwrites DFTI_PACKED_FORMAT
		status = DftiSetValue(desc->hand_back, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
		if (status != 0) { return status; }

		// set configuration: out-of-place
		status = DftiSetValue(desc->hand_back, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
		if (status != 0) { return status; }

		// set complex input strides
		MKL_LONG strides[3];
		strides[0] = 0;
		strides[1] = Nx/2+1;
		strides[2] = 1;
		status = DftiSetValue(desc->hand_back, DFTI_INPUT_STRIDES, strides);
		if (status != 0) { return status; }

		// set real output strides
		strides[0] = 0;
		strides[1] = Nx;
		strides[2] = 1;
		status = DftiSetValue(desc->hand_back, DFTI_OUTPUT_STRIDES, strides);
		if (status != 0) { return status; }

		// commit the descriptor
		status = DftiCommitDescriptor(desc->hand_back);
		if (status != 0) { return status; }
	}

	// normalization factor
	desc->nfac = 1.0 / (Nx * Ny);

	return status;
}


//________________________________________________________________________________________________________________________
///
/// \brief Delete FFT descriptors
///
void DeleteLatticeFourierDescriptor(lattice_fourier_desc_t *desc)
{
	DftiFreeDescriptor(&desc->hand_back);
	DftiFreeDescriptor(&desc->hand_forw);

	desc->Ny = 0;
	desc->Nx = 0;
}


#endif
