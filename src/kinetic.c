#include "kinetic.h"
#include "linalg.h"
#include "util.h"


//________________________________________________________________________________________________________________________
///
/// \brief Allocate and calculate matrix exponential of the kinetic hopping matrix
/// for a rectangular lattice geometry with periodic boundary conditions
///
void RectangularKineticExponential(const sim_params_t *restrict params, kinetic_t *restrict kinetic)
{
	const int Norb = params->Norb;
	const int Nx   = params->Nx;
	const int Ny   = params->Ny;

	const int Ncell = Nx*Ny;
	const int N     = Norb*Ncell;

	kinetic->Norb  = Norb;
	kinetic->Ncell = Ncell;

	int i, j;   // spatial lattice indices
	int o, p;   // orbital indices

	double *T = (double *)algn_calloc(N*N, sizeof(double));

	// set hopping terms in 'T'
	for (o = 0; o < Norb; o++)
	{
		for (p = 0; p < Norb; p++)
		{
			for (j = 0; j < Ny; j++)
			{
				const int j_next = (j < Ny-1 ? j + 1 : 0);

				for (i = 0; i < Nx; i++)
				{
					const int i_next = (i < Nx-1 ? i + 1 : 0);

					// shifted index in x-direction affects 'c' and 'd' cells only
					const int is = (i + (j_next == 0 ? params->pbc_shift : 0)) % Nx;
					const int is_next = (is < Nx-1 ? is + 1 : 0);

					// cell indices
					//  c    d
					//  |\  /
					//  | \/
					//  | /\
					//  |/  \
					//  a----b
					const int a = i       + j     *Nx;
					const int b = i_next  + j     *Nx;
					const int c = is      + j_next*Nx;
					const int d = is_next + j_next*Nx;

					// here, -= is used to avoid overwriting some bonds that are repeated due to
					// periodic boundary conditions, when Nx and/or Ny == 2
					T[(a + o*Ncell)   + (a + p*Ncell)*N] -= params->t.aa[o + p * Norb];
					T[(a + o*Ncell)*N + (a + p*Ncell)  ] -= params->t.aa[o + p * Norb];

					T[(a + o*Ncell)   + (b + p*Ncell)*N] -= params->t.ab[o + p * Norb];
					T[(a + o*Ncell)*N + (b + p*Ncell)  ] -= params->t.ab[o + p * Norb];

					T[(a + o*Ncell)   + (c + p*Ncell)*N] -= params->t.ac[o + p * Norb];
					T[(a + o*Ncell)*N + (c + p*Ncell)  ] -= params->t.ac[o + p * Norb];

					T[(a + o*Ncell)   + (d + p*Ncell)*N] -= params->t.ad[o + p * Norb];
					T[(a + o*Ncell)*N + (d + p*Ncell)  ] -= params->t.ad[o + p * Norb];

					T[(b + o*Ncell)   + (c + p*Ncell)*N] -= params->t.bc[o + p * Norb];
					T[(b + o*Ncell)*N + (c + p*Ncell)  ] -= params->t.bc[o + p * Norb];
				}
			}
		}
	}

	// set diagonal entries in 'T' (shift by chemical potential and site energies)
	for (o = 0; o < Norb; o++)
	{
		for (i = 0; i < Ncell; i++)
		{
			const int index = i + o*Ncell;
			T[index + N*index] = -params->mu + params->eps[o];
		}
	}

	// scale by -dt
	cblas_dscal(N*N, -params->dt, T, 1);

	kinetic->expK     = (double *)algn_malloc(N*N * sizeof(double));
	kinetic->inv_expK = (double *)algn_malloc(N*N * sizeof(double));

	// compute matrix exponential of the kinetic energy operator
	MatrixExp(N, T, kinetic->expK);

	// flip the sign of 'T'
	cblas_dscal(N*N, -1.0, T, 1);

	// compute the inverse matrix exponential of the kinetic energy operator
	MatrixExp(N, T, kinetic->inv_expK);

	// clean up
	algn_free(T);
}


//________________________________________________________________________________________________________________________
///
/// \brief Delete matrix exponentials of kinetic energy operator (free memory)
///
void DeleteKineticExponential(kinetic_t *restrict kinetic)
{
	algn_free(kinetic->inv_expK);
	algn_free(kinetic->expK);
}
