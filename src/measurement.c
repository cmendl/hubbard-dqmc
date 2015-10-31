#include "measurement.h"
#include "util.h"
#include "dupio.h"
#include <mkl.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <assert.h>



//________________________________________________________________________________________________________________________
///
/// \brief Construct the lattice coordinate sum map
///
static void ConstructLatticeCoordinateSumMap(const int Nx, const int Ny, int *restrict latt_sum_map)
{
	// total number of lattice sites
	const int N = Nx * Ny;

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

					latt_sum_map[ij + N*kl] = u + Nx*v;
				}
			}
		}
	}
}


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
	ConstructLatticeCoordinateSumMap(Nx, Ny, meas_data->latt_sum_map);

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

			// k = 0 is special
			meas_data->uu_corr[k] += signfac*(k == 0 ? (1.0 - Gu_ii) : (1.0 - Gu_ii)*(1.0 - Gu_jj) - Gu_ij*Gu_ji);
			meas_data->dd_corr[k] += signfac*(k == 0 ? (1.0 - Gd_ii) : (1.0 - Gd_ii)*(1.0 - Gd_jj) - Gd_ij*Gd_ji);
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
void AllocateUnequalTimeMeasurementData(const int Nx, const int Ny, const int L, measurement_data_unequal_time_t *restrict meas_data)
{
	// total number of lattice sites
	const int N = Nx * Ny;
	meas_data->N = N;

	// total number of time steps
	meas_data->L = L;

	// allocate and initialize Green's functions with zeros
	meas_data->Gtau0_u = (double *)MKL_calloc(L*N*N, sizeof(double), MEM_DATA_ALIGN);
	meas_data->G0tau_u = (double *)MKL_calloc(L*N*N, sizeof(double), MEM_DATA_ALIGN);
	meas_data->Geqlt_u = (double *)MKL_calloc(L*N*N, sizeof(double), MEM_DATA_ALIGN);
	meas_data->Gtau0_d = (double *)MKL_calloc(L*N*N, sizeof(double), MEM_DATA_ALIGN);
	meas_data->G0tau_d = (double *)MKL_calloc(L*N*N, sizeof(double), MEM_DATA_ALIGN);
	meas_data->Geqlt_d = (double *)MKL_calloc(L*N*N, sizeof(double), MEM_DATA_ALIGN);

	// construct lattice coordinate sum map
	meas_data->latt_sum_map = (int *)MKL_malloc(N*N * sizeof(int), MEM_DATA_ALIGN);
	ConstructLatticeCoordinateSumMap(Nx, Ny, meas_data->latt_sum_map);

	// spin correlation data for all time differences
	meas_data->zz_corr = (double *)MKL_calloc(L*N, sizeof(double), MEM_DATA_ALIGN);
	meas_data->xx_corr = (double *)MKL_calloc(L*N, sizeof(double), MEM_DATA_ALIGN);

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

	MKL_free(meas_data->xx_corr);
	MKL_free(meas_data->zz_corr);

	MKL_free(meas_data->latt_sum_map);

	MKL_free(meas_data->Geqlt_d);
	MKL_free(meas_data->G0tau_d);
	MKL_free(meas_data->Gtau0_d);
	MKL_free(meas_data->Geqlt_u);
	MKL_free(meas_data->G0tau_u);
	MKL_free(meas_data->Gtau0_u);
}


//________________________________________________________________________________________________________________________
///
/// \brief Accumulate unequal time "measurement" data
///
void AccumulateUnequalTimeMeasurement(const int sign, const double *const *Bu, const double *const *Bd, measurement_data_unequal_time_t *restrict meas_data)
{
	const int N = meas_data->N;
	const int L = meas_data->L;

	// spin-up
	double *curGtau0_u = (double *)MKL_malloc(L*N*N * sizeof(double), MEM_DATA_ALIGN);
	double *curG0tau_u = (double *)MKL_malloc(L*N*N * sizeof(double), MEM_DATA_ALIGN);
	double *curGeqlt_u = (double *)MKL_malloc(L*N*N * sizeof(double), MEM_DATA_ALIGN);

	// spin-down
	double *curGtau0_d = (double *)MKL_malloc(L*N*N * sizeof(double), MEM_DATA_ALIGN);
	double *curG0tau_d = (double *)MKL_malloc(L*N*N * sizeof(double), MEM_DATA_ALIGN);
	double *curGeqlt_d = (double *)MKL_malloc(L*N*N * sizeof(double), MEM_DATA_ALIGN);

	// compute unequal time spin-up and spin-down Green's functions
	#pragma omp parallel sections
	{
		#pragma omp section
		ComputeUnequalTimeGreensFunction(N, L, Bu, meas_data->Hu, curGtau0_u, curG0tau_u, curGeqlt_u, NULL, NULL, NULL);
		#pragma omp section
		ComputeUnequalTimeGreensFunction(N, L, Bd, meas_data->Hd, curGtau0_d, curG0tau_d, curGeqlt_d, NULL, NULL, NULL);
	}

	// accumulate unequal time Green's functions
	cblas_daxpy(L*N*N, sign, curGtau0_u, 1, meas_data->Gtau0_u, 1);
	cblas_daxpy(L*N*N, sign, curG0tau_u, 1, meas_data->G0tau_u, 1);
	cblas_daxpy(L*N*N, sign, curGeqlt_u, 1, meas_data->Geqlt_u, 1);
	cblas_daxpy(L*N*N, sign, curGtau0_d, 1, meas_data->Gtau0_d, 1);
	cblas_daxpy(L*N*N, sign, curG0tau_d, 1, meas_data->G0tau_d, 1);
	cblas_daxpy(L*N*N, sign, curGeqlt_d, 1, meas_data->Geqlt_d, 1);

	// TODO: accumulate spin correlation data

	// clean up
	MKL_free(curGeqlt_d);
	MKL_free(curG0tau_d);
	MKL_free(curGtau0_d);
	MKL_free(curGeqlt_u);
	MKL_free(curG0tau_u);
	MKL_free(curGtau0_u);

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

	// divide Green's functions by sign
	cblas_dscal(L*N*N, nfac, meas_data->Gtau0_u, 1);
	cblas_dscal(L*N*N, nfac, meas_data->G0tau_u, 1);
	cblas_dscal(L*N*N, nfac, meas_data->Geqlt_u, 1);
	cblas_dscal(L*N*N, nfac, meas_data->Gtau0_d, 1);
	cblas_dscal(L*N*N, nfac, meas_data->G0tau_d, 1);
	cblas_dscal(L*N*N, nfac, meas_data->Geqlt_d, 1);

	// calculate average sign
	meas_data->sign /= meas_data->nsampl;

	// set sample counter to 1
	meas_data->nsampl = 1;
}
