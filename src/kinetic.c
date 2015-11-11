#include "kinetic.h"
#include "linalg.h"
#include <mkl.h>


//________________________________________________________________________________________________________________________
///
/// \brief Allocate and calculate matrix exponential of the kinetic nearest and next-nearest
/// neighbor hopping matrix for a rectangular lattice geometry with periodic boundary conditions.
///
void RectangularKineticExponential(const sim_params_t *restrict params, kinetic_t *restrict kinetic)
{
	const int Norb = params->Norb;
	const int Nx = params->Nx;
	const int Ny = params->Ny;
	const int N = Norb*Nx*Ny;

	kinetic->Norb = Norb;
	kinetic->Ncell = Nx*Ny;
	kinetic->N = N;

	int i, j; // spatial indices
	int o, p; // these are orbital numbers

	double *T = (double *)MKL_calloc(N*N, sizeof(double), MEM_DATA_ALIGN);

	// set hopping terms in 'T'
	for (j = 0; j < Ny; j++)
	{
		const int j_next = (j < Ny-1 ? j + 1 : 0);

		for (i = 0; i < Nx; i++)
		{
			const int i_next = (i < Nx-1 ? i + 1 : 0);

			for (o = 0; o < Norb; o++)
			{
				for (p = 0; p < Norb; p++)
				{
					T[(i + j*Nx + o*Nx*Ny) + N*(i      + j     *Nx + p*Nx*Ny)] = params->t.aa[o + p * Norb];
					T[(i + j*Nx + o*Nx*Ny) + N*(i_next + j     *Nx + p*Nx*Ny)] = params->t.ab[o + p * Norb];
					T[(i + j*Nx + o*Nx*Ny) + N*(i      + j_next*Nx + p*Nx*Ny)] = params->t.ac[o + p * Norb];
					T[(i + j*Nx + o*Nx*Ny) + N*(i_next + j_next*Nx + p*Nx*Ny)] = params->t.ad[o + p * Norb];
					T[(i_next + j*Nx + o*Nx*Ny) + N*(i + j_next*Nx + p*Nx*Ny)] = params->t.bc[o + p * Norb];

					T[(i      + j     *Nx + p*Nx*Ny) + (i + j*Nx + o*Nx*Ny)*N] = params->t.aa[o + p * Norb];
					T[(i_next + j     *Nx + p*Nx*Ny) + (i + j*Nx + o*Nx*Ny)*N] = params->t.ab[o + p * Norb];
					T[(i      + j_next*Nx + p*Nx*Ny) + (i + j*Nx + o*Nx*Ny)*N] = params->t.ac[o + p * Norb];
					T[(i_next + j_next*Nx + p*Nx*Ny) + (i + j*Nx + o*Nx*Ny)*N] = params->t.ad[o + p * Norb];
					T[(i + j_next*Nx + p*Nx*Ny) + (i_next + j*Nx + p*Nx*Ny)*N] = params->t.bc[o + p * Norb];
				}
			}
		}
	}

	// set diagonal entries in 'T' (shift by chemical potential and site energies)
	for (o = 0; o < Norb; o++)
	{
		for (i = 0; i < Nx*Ny; i++)
		{
			const int index = i + o*Nx*Ny;
			T[index + N * index] = params->mu + params->eps[o];
		}
	}

	//scale by dt
	cblas_dscal(N*N, params->dt, T, 1);

	kinetic->expK     = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);
	kinetic->inv_expK = (double *)MKL_malloc(N*N * sizeof(double), MEM_DATA_ALIGN);

	// compute matrix exponential of the kinetic energy operator
	MatrixExp(N, T, kinetic->expK);

	// flip the sign of 'T'
	cblas_dscal(N*N, -1.0, T, 1);

	// compute the inverse matrix exponential of the kinetic energy operator
	MatrixExp(N, T, kinetic->inv_expK);

	// clean up
	MKL_free(T);
}


//________________________________________________________________________________________________________________________
///
/// \brief Delete matrix exponentials of kinetic energy operator (free memory)
///
void DeleteKineticExponential(kinetic_t *restrict kinetic)
{
	MKL_free(kinetic->inv_expK);
	MKL_free(kinetic->expK);
}
