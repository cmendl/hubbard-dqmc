#include "measurement.h"
#include "linalg.h"
#include "dupio.h"
#include <mkl.h>
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

	int i;
	for (i = 0; i < 2; i++)
	{
		meas_data->density_u[i] = 0;
		meas_data->density_d[i] = 0;
		meas_data->doubleocc[i] = 0;
	}

	// construct lattice coordinate sum map
	meas_data->latt_sum_map = (int *)MKL_malloc(N*N * sizeof(int), MEM_DATA_ALIGN);
	int j;
	for (j = 0; j < Ny; j++)
	{
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
void AccumulateEqualTimeMeasurement(const greens_func_t *restrict Gu, const greens_func_t *restrict Gd, measurement_data_t *meas_data)
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

	// mean value
	meas_data->density_u[0] += sign*nu;
	meas_data->density_d[0] += sign*nd;
	meas_data->doubleocc[0] += sign*oc;
	// squares
	meas_data->density_u[1] += square(nu);
	meas_data->density_d[1] += square(nd);
	meas_data->doubleocc[1] += square(oc);

	// density and spin correlations
	for (i = 0; i < N; i++)
	{
		const double Gu_ii = Gu->mat[i + N*i];
		const double Gd_ii = Gd->mat[i + N*i];

		for (k = 0; k < N; k++)
		{
			j = meas_data->latt_sum_map[k + N*i];

			const double Gu_jj = Gu->mat[j + N*j];
			const double Gd_jj = Gd->mat[j + N*j];
			const double Gu_ij = Gu->mat[i + N*j];
			const double Gd_ij = Gd->mat[i + N*j];
			const double Gu_ji = Gu->mat[j + N*i];
			const double Gd_ji = Gd->mat[j + N*i];

			meas_data->uu_corr[k] += sign*((1.0 - Gu_ii)*(1.0 - Gu_jj) - Gu_ij*Gu_ji);
			meas_data->dd_corr[k] += sign*((1.0 - Gd_ii)*(1.0 - Gd_jj) - Gd_ij*Gd_ji);
			meas_data->ud_corr[k] += sign*((1.0 - Gu_ii)*(1.0 - Gd_jj) + (1.0 - Gd_ii)*(1.0 - Gu_jj));

			meas_data->zz_corr[k] += sign*((Gu_ii - Gd_ii)*(Gu_jj - Gd_jj) - (Gu_ij*Gu_ji + Gd_ij*Gd_ji));
			meas_data->xx_corr[k] += sign*(                                - (Gu_ij*Gd_ji + Gd_ij*Gu_ji));
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

	// replace (close to) zero sign by 1 to avoid division by zero
	if (fabs(meas_data->sign) < 0.5)
	{
		duprintf("Warning: accumulated sign '%g' is (close to) zero, replacing it by 1 (%i samples).\n", meas_data->sign, meas_data->nsampl);
		meas_data->sign = 1.0;
	}

	// normalization factor
	const double nfac = 1.0 / meas_data->sign;

	int i;
	for (i = 0; i < 2; i++)
	{
		meas_data->density_u[i] *= nfac;
		meas_data->density_d[i] *= nfac;
		meas_data->doubleocc[i] *= nfac;
	}

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
