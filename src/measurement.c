#include "measurement.h"
#include "linalg.h"
#include "util.h"
#include "dupio.h"
#include <mkl.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <assert.h>


//________________________________________________________________________________________________________________________
///
/// \brief Allocate and initialize measurement data structure
///
void AllocateMeasurementData(const int Nx, const int Ny, measurement_data_t *restrict meas_data)
{
	// total number of lattice sites
	const int N = Nx * Ny;
	meas_data->N = N;

	meas_data->density_u = 0;
	meas_data->density_d = 0;
	meas_data->doubleocc = 0;

	// construct lattice coordinate sum map
	meas_data->latt_sum_map = (int *)MKL_malloc(N*N * sizeof(int), MEM_DATA_ALIGN);
	int j;
	for (j = 0; j < Ny; j++)
	{
		int i;
		for (i = 0; i < Nx; i++)
		{
			// index of (i,j) lattice site
			const int ij = i + Nx*j;

			int k, l;
			for (l = 0; l < Ny; l++)
			{
				for (k = 0; k < Nx; k++)
				{
					// index of (k,l) lattice site
					const int kl = k + Nx*l;

					// coordinate sum
					const int u = (i + k) % Nx;
					const int v = (j + l) % Ny;

					meas_data->latt_sum_map[ij + N*kl] = u + Nx*v;
				}
			}
		}
	}

	// density correlation data
	meas_data->uu_corr = (double *)MKL_calloc(N, sizeof(double), MEM_DATA_ALIGN);
	meas_data->dd_corr = (double *)MKL_calloc(N, sizeof(double), MEM_DATA_ALIGN);
	meas_data->ud_corr = (double *)MKL_calloc(N, sizeof(double), MEM_DATA_ALIGN);

	// spin correlation data
	meas_data->zz_corr = (double *)MKL_calloc(N, sizeof(double), MEM_DATA_ALIGN);
	meas_data->xx_corr = (double *)MKL_calloc(N, sizeof(double), MEM_DATA_ALIGN);

	meas_data->sign = 0;

	// no samples collected so far
	meas_data->nsampl = 0;
}


//________________________________________________________________________________________________________________________
///
/// \brief Delete measurement data (free memory)
///
void DeleteMeasurementData(measurement_data_t *restrict meas_data)
{
	MKL_free(meas_data->xx_corr);
	MKL_free(meas_data->zz_corr);
	MKL_free(meas_data->ud_corr);
	MKL_free(meas_data->dd_corr);
	MKL_free(meas_data->uu_corr);
	MKL_free(meas_data->latt_sum_map);
}


//________________________________________________________________________________________________________________________
///
/// \brief Accumulate equal time "measurement" data
///
void AccumulateMeasurement(const greens_func_t *restrict Gu, const greens_func_t *restrict Gd, measurement_data_t *restrict meas_data)
{
	int i, j, k;

	// total number of lattice sites
	const int N = meas_data->N;

	// product of the determinant signs of the Green's function matrices
	const double sign = (double)(Gu->sgndet * Gd->sgndet);
	assert(sign != 0);

	double nu = 0;
	double nd = 0;
	double oc = 0;
	for (i = 0; i < N; i++)
	{
		nu += (1 - Gu->mat[i + i*N]);
		nd += (1 - Gd->mat[i + i*N]);
		oc += (1 - Gu->mat[i + i*N])*(1 - Gd->mat[i + i*N]);
	}
	// normalization
	const double nfac = 1.0 / N;
	nu *= nfac;
	nd *= nfac;
	oc *= nfac;

	// mean value of density and double occupancy
	meas_data->density_u += sign*nu;
	meas_data->density_d += sign*nd;
	meas_data->doubleocc += sign*oc;

	// density and spin correlations

	// sign and normalization factor
	const double signfac = sign / N;

	for (i = 0; i < N; i++)
	{
		const double Gu_ii = Gu->mat[i + N*i];
		const double Gd_ii = Gd->mat[i + N*i];

		meas_data->zz_corr[0] += signfac*(Gu_ii + Gd_ii);
		meas_data->xx_corr[0] += signfac*(Gu_ii + Gd_ii);

		for (k = 0; k < N; k++)
		{
			j = meas_data->latt_sum_map[k + N*i];

			const double Gu_jj = Gu->mat[j + N*j];
			const double Gd_jj = Gd->mat[j + N*j];
			const double Gu_ij = Gu->mat[i + N*j];
			const double Gd_ij = Gd->mat[i + N*j];
			const double Gu_ji = Gu->mat[j + N*i];
			const double Gd_ji = Gd->mat[j + N*i];

			meas_data->uu_corr[k] += signfac*((1.0 - Gu_ii)*(1.0 - Gu_jj) - Gu_ij*Gu_ji);
			meas_data->dd_corr[k] += signfac*((1.0 - Gd_ii)*(1.0 - Gd_jj) - Gd_ij*Gd_ji);
			meas_data->ud_corr[k] += signfac*((1.0 - Gu_ii)*(1.0 - Gd_jj) + (1.0 - Gd_ii)*(1.0 - Gu_jj));

			meas_data->zz_corr[k] += signfac*((Gu_ii - Gd_ii)*(Gu_jj - Gd_jj) - (Gu_ij*Gu_ji + Gd_ij*Gd_ji));
			meas_data->xx_corr[k] += signfac*(                                - (Gu_ij*Gd_ji + Gd_ij*Gu_ji));
		}
	}

	// add current sign
	meas_data->sign += sign;

	// increment sample counter
	meas_data->nsampl++;
}


//________________________________________________________________________________________________________________________
///
/// \brief Normalize measurement data (divide by accumulated sign)
///
void NormalizeMeasurementData(measurement_data_t *meas_data)
{
	// total number of lattice sites
	const int N = meas_data->N;

	// normalization factor; sign must be non-zero
	const double nfac = 1.0 / meas_data->sign;

	meas_data->density_u *= nfac;
	meas_data->density_d *= nfac;
	meas_data->doubleocc *= nfac;

	int i;
	for (i = 0; i < N; i++)
	{
		meas_data->uu_corr[i] *= nfac;
		meas_data->dd_corr[i] *= nfac;
		meas_data->ud_corr[i] *= nfac;

		meas_data->zz_corr[i] *= nfac;
		meas_data->xx_corr[i] *= nfac;
	}

	// calculate average sign
	meas_data->sign /= meas_data->nsampl;

	// set sample counter to 1
	meas_data->nsampl = 1;
}


//________________________________________________________________________________________________________________________
///
/// \brief Allocate and initialize unequal time measurement data structure
///
void AllocateUnequalTimeMeasurementData(const int N, const int L, measurement_data_unequal_time_t *restrict meas_data)
{
	meas_data->N = N;
	meas_data->L = L;

	// initialize with zeros
	meas_data->Gu_ut = (double *)MKL_calloc(L*N*N, sizeof(double), MEM_DATA_ALIGN);
	meas_data->Gd_ut = (double *)MKL_calloc(L*N*N, sizeof(double), MEM_DATA_ALIGN);

	meas_data->Hu = (double *)MKL_malloc(N*N*L*L * sizeof(double), MEM_DATA_ALIGN);
	meas_data->Hd = (double *)MKL_malloc(N*N*L*L * sizeof(double), MEM_DATA_ALIGN);

	meas_data->sign = 0;

	// no samples collected so far
	meas_data->nsampl = 0;
}


//________________________________________________________________________________________________________________________
///
/// \brief Delete unequal time measurement data (free memory)
///
void DeleteUnequalTimeMeasurementData(measurement_data_unequal_time_t *restrict meas_data)
{
	MKL_free(meas_data->Hd);
	MKL_free(meas_data->Hu);

	MKL_free(meas_data->Gd_ut);
	MKL_free(meas_data->Gu_ut);
}


//________________________________________________________________________________________________________________________
///
/// \brief Accumulate unequal time Green's function
///
/// Reference:
///   - R. Blankenbecler, D. J. Scalapino, R. L. Sugar\n
///     Monte Carlo calculations of coupled boson-fermion systems. I\n
///     Phys. Rev. D 24, 2278 (1981)
///
static void AccumulateUnequalTimeGreensFunction(const int sign, const int N, const int L, const double *const *B, double *restrict H, double *restrict G)
{
	int i, j, k;

	// set H to the identity matrix
	memset(H, 0, N*N*L*L * sizeof(double));
	for (i = 0; i < L*N; i++)
	{
		H[i + i*L*N] = 1;
	}

	// copy B matrices as blocks into H
	for (i = 0; i < L; i++)
	{
		if (i == 0)
		{
			// upper right block
			mkl_domatcopy('C', 'N', N, N, 1.0, B[i], N, &H[L*N*(L - 1)*N], L*N);
		}
		else
		{
			// lower off-diagonal blocks; note the factor (-1)
			mkl_domatcopy('C', 'N', N, N, -1.0, B[i], N, &H[(i + L*N*(i - 1))*N], L*N);
		}
	}

	// compute the inverse of the block L-cyclic H matrix
	BlockCyclicInverse(N, L, H);

	// sign and normalization factor
	const double signfac = (double)sign / L;

	// accumulate cyclically shifted column blocks of H^{-1} into G
	for (j = 0; j < L; j++)
	{
		for (i = 0; i < L; i++)
		{
			const int ij = (i + j) % L;

			for (k = 0; k < N; k++)
			{
				cblas_daxpy(N, signfac, &H[(ij + L*(k + N*j))*N], 1, &G[(i + k*L)*N], 1);
			}
		}
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief Accumulate unequal time "measurement" data
///
void AccumulateUnequalTimeMeasurement(const int sign, const double *const *Bu, const double *const *Bd, measurement_data_unequal_time_t *restrict meas_data)
{
	#pragma omp parallel sections
	{
		#pragma omp section
		AccumulateUnequalTimeGreensFunction(sign, meas_data->N, meas_data->L, Bu, meas_data->Hu, meas_data->Gu_ut);
		#pragma omp section
		AccumulateUnequalTimeGreensFunction(sign, meas_data->N, meas_data->L, Bd, meas_data->Hd, meas_data->Gd_ut);
	}

	// add current sign
	meas_data->sign += (double)sign;

	// increment sample counter
	meas_data->nsampl++;
}


//________________________________________________________________________________________________________________________
///
/// \brief Normalize unequal time measurement data (divide by accumulated sign)
///
void NormalizeUnequalTimeMeasurementData(measurement_data_unequal_time_t *meas_data)
{
	// total number of lattice sites
	const int N = meas_data->N;

	// total number of time steps
	const int L = meas_data->L;

	// normalization factor; sign must be non-zero
	const double nfac = 1.0 / meas_data->sign;

	cblas_dscal(L*N*N, nfac, meas_data->Gu_ut, 1);
	cblas_dscal(L*N*N, nfac, meas_data->Gd_ut, 1);

	// calculate average sign
	meas_data->sign /= meas_data->nsampl;

	// set sample counter to 1
	meas_data->nsampl = 1;
}
