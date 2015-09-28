#include "kinetic.h"
#include "linalg.h"
#include <mkl.h>


//________________________________________________________________________________________________________________________
///
/// \brief Allocate and calculate matrix exponential of the kinetic nearest neighbor hopping matrix
///
void NearestNeighborKineticExponential(const int Nx, const int Ny, const double mu, const double dt, kinetic_t *restrict kinetic)
{
	// total number of lattice sites
	const int N = Nx * Ny;
	kinetic->N = N;

	int i, j;

	double *T = (double *)MKL_calloc(N*N, sizeof(double), MEM_DATA_ALIGN);

	// set diagonal entries in 'T' (shift by chemical potential)
	const double dt_mu = dt * mu;
	for (i = 0; i < Nx; i++)
	{
		T[i + N*i] = dt_mu;
	}

	// set hopping terms in 'T'
	for (j = 0; j < Ny; j++)
	{
		const int j_next = (j < Ny-1 ? j + 1 : 0   );
		const int j_prev = (j > 0    ? j - 1 : Ny-1);

		for (i = 0; i < Nx; i++)
		{
			const int i_next = (i < Nx-1 ? i + 1 : 0   );
			const int i_prev = (i > 0    ? i - 1 : Nx-1);

			T[(i_next + j*Nx) + N*(i      + j*Nx)] = dt;
			T[(i      + j*Nx) + N*(i_next + j*Nx)] = dt;

			T[(i_prev + j*Nx) + N*(i      + j*Nx)] = dt;
			T[(i      + j*Nx) + N*(i_prev + j*Nx)] = dt;

			T[(i + j_next*Nx) + N*(i + j     *Nx)] = dt;
			T[(i + j     *Nx) + N*(i + j_next*Nx)] = dt;

			T[(i + j_prev*Nx) + N*(i + j     *Nx)] = dt;
			T[(i + j     *Nx) + N*(i + j_prev*Nx)] = dt;
		}
	}

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
