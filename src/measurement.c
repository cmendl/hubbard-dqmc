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
	const int Ncell = Nx * Ny;

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

					latt_sum_map[ij + Ncell*kl] = u + Nx*v;
				}
			}
		}
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief Construct lattice nearest neighbor maps
///
static void ConstructLatticeNearestNeighborMap(const int Nx, const int Ny, int *restrict latt_xp1_map, int *restrict latt_xm1_map, int *restrict latt_yp1_map, int *restrict latt_ym1_map)
{
	int j;
	for (j = 0; j < Ny; j++)
	{
		const int j_next = (j < Ny-1 ? j + 1 : 0   );
		const int j_prev = (j > 0    ? j - 1 : Ny-1);

		int i;
		for (i = 0; i < Nx; i++)
		{
			const int i_next = (i < Nx-1 ? i + 1 : 0   );
			const int i_prev = (i > 0    ? i - 1 : Nx-1);

			// index of (i,j) lattice site
			const int ij = i + Nx*j;

			latt_xp1_map[ij] = i_next + Nx*j;			// index of (i+1,j) lattice site
			latt_xm1_map[ij] = i_prev + Nx*j;			// index of (i-1,j) lattice site
			latt_yp1_map[ij] = i      + Nx*j_next;		// index of (i,j+1) lattice site
			latt_ym1_map[ij] = i      + Nx*j_prev;		// index of (i,j-1) lattice site
		}
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief Allocate and initialize measurement data structure
///
void AllocateMeasurementData(const int Norb, const int Nx, const int Ny, measurement_data_t *restrict meas_data)
{
	const int Ncell = Nx * Ny;

	meas_data->Norb  = Norb;
	meas_data->Ncell = Ncell;

	// no samples collected so far
	meas_data->nsampl = 0;

	// construct lattice coordinate sum map
	meas_data->latt_sum_map = (int *)MKL_malloc(Ncell*Ncell * sizeof(int), MEM_DATA_ALIGN);
	ConstructLatticeCoordinateSumMap(Nx, Ny, meas_data->latt_sum_map);

	meas_data->sign = 0;

	meas_data->density_u = (double *)MKL_calloc(Norb, sizeof(double), MEM_DATA_ALIGN);
	meas_data->density_d = (double *)MKL_calloc(Norb, sizeof(double), MEM_DATA_ALIGN);
	meas_data->doubleocc = (double *)MKL_calloc(Norb, sizeof(double), MEM_DATA_ALIGN);

	// density correlation data
	// first index for spatial offset, second and third indices for orbital index
	meas_data->uu_corr = (double *)MKL_calloc(Ncell*Norb*Norb, sizeof(double), MEM_DATA_ALIGN);
	meas_data->dd_corr = (double *)MKL_calloc(Ncell*Norb*Norb, sizeof(double), MEM_DATA_ALIGN);
	meas_data->ud_corr = (double *)MKL_calloc(Ncell*Norb*Norb, sizeof(double), MEM_DATA_ALIGN);
	meas_data->ff_corr = (double *)MKL_calloc(Ncell*Norb*Norb, sizeof(double), MEM_DATA_ALIGN);

	// spin correlation data
	meas_data->zz_corr = (double *)MKL_calloc(Ncell*Norb*Norb, sizeof(double), MEM_DATA_ALIGN);
	meas_data->xx_corr = (double *)MKL_calloc(Ncell*Norb*Norb, sizeof(double), MEM_DATA_ALIGN);
}


//________________________________________________________________________________________________________________________
///
/// \brief Delete measurement data (free memory)
///
void DeleteMeasurementData(measurement_data_t *restrict meas_data)
{
	MKL_free(meas_data->xx_corr);
	MKL_free(meas_data->zz_corr);
	MKL_free(meas_data->ff_corr);
	MKL_free(meas_data->ud_corr);
	MKL_free(meas_data->dd_corr);
	MKL_free(meas_data->uu_corr);
	MKL_free(meas_data->doubleocc);
	MKL_free(meas_data->density_d);
	MKL_free(meas_data->density_u);
	MKL_free(meas_data->latt_sum_map);
}


//________________________________________________________________________________________________________________________
///
/// \brief Accumulate equal time "measurement" data
///
void AccumulateMeasurement(const greens_func_t *restrict Gu, const greens_func_t *restrict Gd, measurement_data_t *restrict meas_data)
{
	int i, k;	// spatial indices
	int o, p;	// orbital indices

	// lattice dimensions
	const int Norb  = meas_data->Norb;
	const int Ncell = meas_data->Ncell;
	const int N = Norb * Ncell;

	// product of the determinant signs of the Green's function matrices
	const double sign = (double)(Gu->sgndet * Gd->sgndet);
	assert(sign != 0);

	// sign and normalization factor
	const double signfac = sign / Ncell;

	for (o = 0; o < Norb; o++)
	{
		double nu = 0;
		double nd = 0;
		double oc = 0;
		for (i = 0; i < Ncell; i++)
		{
			const int io = i + o * Ncell;
			nu += (1 - Gu->mat[io + io*N]);
			nd += (1 - Gd->mat[io + io*N]);
			oc += (1 - Gu->mat[io + io*N])*(1 - Gd->mat[io + io*N]);
		}
		// mean value of density and double occupancy for current orbital
		meas_data->density_u[o] += signfac*nu;
		meas_data->density_d[o] += signfac*nd;
		meas_data->doubleocc[o] += signfac*oc;
	}

	// density and spin correlations

	for (o = 0; o < Norb; o++)
	{
		for (p = 0; p < Norb; p++)
		{
			// starting index for correlation accumulators
			const int offset = (o + p*Norb) * Ncell;

			for (i = 0; i < Ncell; i++)
			{
				const int io = i + o * Ncell;
				const double Gu_ii = Gu->mat[io + N*io];
				const double Gd_ii = Gd->mat[io + N*io];

				for (k = 0; k < Ncell; k++)
				{
					const int jp = meas_data->latt_sum_map[k + Ncell*i] + p * Ncell;

					const double Gu_jj = Gu->mat[jp + N*jp];
					const double Gd_jj = Gd->mat[jp + N*jp];
					const double Gu_ij = Gu->mat[io + N*jp];
					const double Gd_ij = Gd->mat[io + N*jp];
					const double Gu_ji = Gu->mat[jp + N*io];
					const double Gd_ji = Gd->mat[jp + N*io];

					if (k == 0 && o == p)	// special case: same lattice site and orbital
					{
						meas_data->uu_corr[k + offset] += signfac*(1 - Gu_ii);
						meas_data->dd_corr[k + offset] += signfac*(1 - Gd_ii);
						meas_data->ff_corr[k + offset] += signfac*(Gu_ii*(1 - Gu_ii) + Gd_ii*(1 - Gd_ii));
						meas_data->zz_corr[k + offset] += signfac*(Gu_ii - 2*Gu_ii*Gd_ii + Gd_ii);
						meas_data->xx_corr[k + offset] += signfac*(Gu_ii - 2*Gu_ii*Gd_ii + Gd_ii);
					}
					else
					{
						meas_data->uu_corr[k + offset] += signfac*((1 - Gu_ii)*(1 - Gu_jj) - Gu_ij*Gu_ji);
						meas_data->dd_corr[k + offset] += signfac*((1 - Gd_ii)*(1 - Gd_jj) - Gd_ij*Gd_ji);
						meas_data->ff_corr[k + offset] += signfac*(                                - (Gu_ij*Gu_ji + Gd_ij*Gd_ji));
						meas_data->zz_corr[k + offset] += signfac*((Gu_ii - Gd_ii)*(Gu_jj - Gd_jj) - (Gu_ij*Gu_ji + Gd_ij*Gd_ji));
						meas_data->xx_corr[k + offset] += signfac*(                                - (Gu_ij*Gd_ji + Gd_ij*Gu_ji));
					}

					// independent of special case
					meas_data->ud_corr[k + offset] += signfac*((1 - Gu_ii)*(1 - Gd_jj) + (1 - Gd_ii)*(1 - Gu_jj));
				}
			}
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
	// total number of orbitals and cells
	const int Norb = meas_data->Norb;
	const int Ncell = meas_data->Ncell;

	// normalization factor; sign must be non-zero
	const double nfac = 1.0 / meas_data->sign;

	int i;
	for (i = 0; i < Norb; i++)
	{
		meas_data->density_u[i] *= nfac;
		meas_data->density_d[i] *= nfac;
		meas_data->doubleocc[i] *= nfac;
	}

	for (i = 0; i < Ncell*Norb*Norb; i++)
	{
		meas_data->uu_corr[i] *= nfac;
		meas_data->dd_corr[i] *= nfac;
		meas_data->ud_corr[i] *= nfac;
		meas_data->ff_corr[i] *= nfac;
		meas_data->zz_corr[i] *= nfac;
		meas_data->xx_corr[i] *= nfac;
	}

	// calculate average sign
	meas_data->sign /= meas_data->nsampl;
}


//________________________________________________________________________________________________________________________
///
/// \brief Print a summary of basic measurements
///
void PrintMeasurementDataSummary(const measurement_data_t *meas_data)
{
	duprintf("_______________________________________________________________________________\n");
	duprintf("Summary of simulation results\n\n");
	duprintf("                    average sign: %g\n", meas_data->sign);
	double total_density = 0;
	int o;
	for (o = 0; o < meas_data->Norb; o++)
	{
		total_density += meas_data->density_u[o] + meas_data->density_d[o];
	}
	duprintf("           average total density: %g\n", total_density);
	for (o = 0; o < meas_data->Norb; o++)
	{
		duprintf("\nResults for orbital %d\n", o);
		duprintf("           average total density: %g\n", meas_data->density_u[o] + meas_data->density_d[o]);
		duprintf("         average spin-up density: %g\n", meas_data->density_u[o]);
		duprintf("       average spin-down density: %g\n", meas_data->density_d[o]);
		duprintf("        average double occupancy: %g\n", meas_data->doubleocc[o]);
		duprintf("            average local moment: %g\n", meas_data->density_u[o] + meas_data->density_d[o] - 2.0*meas_data->doubleocc[o]);
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief Load equal time measurement data structure from disk
///
void LoadMeasurementData(const char *fnbase, measurement_data_t *meas_data)
{
	const int Norb  = meas_data->Norb;
	const int Ncell = meas_data->Ncell;

	char path[1024];
	sprintf(path, "%s_sign.dat",      fnbase); ReadData(path, (void *)&meas_data->sign,   sizeof(double), 1);
	sprintf(path, "%s_nsampl.dat",    fnbase); ReadData(path, (void *)&meas_data->nsampl, sizeof(int),    1);
	sprintf(path, "%s_density_u.dat", fnbase); ReadData(path, meas_data->density_u, sizeof(double), Norb);
	sprintf(path, "%s_density_d.dat", fnbase); ReadData(path, meas_data->density_d, sizeof(double), Norb);
	sprintf(path, "%s_doubleocc.dat", fnbase); ReadData(path, meas_data->doubleocc, sizeof(double), Norb);
	sprintf(path, "%s_uu_corr.dat",   fnbase); ReadData(path, meas_data->uu_corr,   sizeof(double), Ncell*Norb*Norb);
	sprintf(path, "%s_dd_corr.dat",   fnbase); ReadData(path, meas_data->dd_corr,   sizeof(double), Ncell*Norb*Norb);
	sprintf(path, "%s_ud_corr.dat",   fnbase); ReadData(path, meas_data->ud_corr,   sizeof(double), Ncell*Norb*Norb);
	sprintf(path, "%s_ff_corr.dat",   fnbase); ReadData(path, meas_data->ff_corr,   sizeof(double), Ncell*Norb*Norb);
	sprintf(path, "%s_zz_corr.dat",   fnbase); ReadData(path, meas_data->zz_corr,   sizeof(double), Ncell*Norb*Norb);
	sprintf(path, "%s_xx_corr.dat",   fnbase); ReadData(path, meas_data->xx_corr,   sizeof(double), Ncell*Norb*Norb);
}


//________________________________________________________________________________________________________________________
///
/// \brief Save equal time measurement data structure to disk
///
void SaveMeasurementData(const char *fnbase, const measurement_data_t *meas_data)
{
	const int Norb  = meas_data->Norb;
	const int Ncell = meas_data->Ncell;

	char path[1024];
	sprintf(path, "%s_sign.dat",      fnbase); WriteData(path, &meas_data->sign,     sizeof(double), 1, false);
	sprintf(path, "%s_nsampl.dat",    fnbase); WriteData(path, &meas_data->nsampl,   sizeof(int),    1, false);
	sprintf(path, "%s_density_u.dat", fnbase); WriteData(path, meas_data->density_u, sizeof(double), Norb, false);
	sprintf(path, "%s_density_d.dat", fnbase); WriteData(path, meas_data->density_d, sizeof(double), Norb, false);
	sprintf(path, "%s_doubleocc.dat", fnbase); WriteData(path, meas_data->doubleocc, sizeof(double), Norb, false);
	sprintf(path, "%s_uu_corr.dat",   fnbase); WriteData(path, meas_data->uu_corr,   sizeof(double), Ncell*Norb*Norb, false);
	sprintf(path, "%s_dd_corr.dat",   fnbase); WriteData(path, meas_data->dd_corr,   sizeof(double), Ncell*Norb*Norb, false);
	sprintf(path, "%s_ud_corr.dat",   fnbase); WriteData(path, meas_data->ud_corr,   sizeof(double), Ncell*Norb*Norb, false);
	sprintf(path, "%s_ff_corr.dat",   fnbase); WriteData(path, meas_data->ff_corr,   sizeof(double), Ncell*Norb*Norb, false);
	sprintf(path, "%s_zz_corr.dat",   fnbase); WriteData(path, meas_data->zz_corr,   sizeof(double), Ncell*Norb*Norb, false);
	sprintf(path, "%s_xx_corr.dat",   fnbase); WriteData(path, meas_data->xx_corr,   sizeof(double), Ncell*Norb*Norb, false);
}


//________________________________________________________________________________________________________________________
///
/// \brief Allocate and initialize unequal time measurement data structure
///
int AllocateUnequalTimeMeasurementData(const int Norb, const int Nx, const int Ny, const int L, measurement_data_unequal_time_t *restrict meas_data)
{
	// lattice dimensions
	const int Ncell = Nx * Ny;
	const int N     = Norb * Ncell;

	meas_data->Norb  = Norb;
	meas_data->Ncell = Ncell;

	// total number of time steps
	meas_data->L = L;

	// allocate temporary H matrices
	meas_data->Hu = (double *)MKL_malloc(L*L*N*N * sizeof(double), MEM_DATA_ALIGN); if (meas_data->Hu == NULL) { return -1; }
	meas_data->Hd = (double *)MKL_malloc(L*L*N*N * sizeof(double), MEM_DATA_ALIGN); if (meas_data->Hd == NULL) { return -1; }

	// allocate and initialize Green's functions with zeros
	meas_data->Gtau0_u = (double *)MKL_calloc(L*N*N, sizeof(double), MEM_DATA_ALIGN); if (meas_data->Gtau0_u == NULL) { return -1; }
	meas_data->G0tau_u = (double *)MKL_calloc(L*N*N, sizeof(double), MEM_DATA_ALIGN); if (meas_data->G0tau_u == NULL) { return -1; }
	meas_data->Geqlt_u = (double *)MKL_calloc(L*N*N, sizeof(double), MEM_DATA_ALIGN); if (meas_data->Geqlt_u == NULL) { return -1; }
	meas_data->Gtau0_d = (double *)MKL_calloc(L*N*N, sizeof(double), MEM_DATA_ALIGN); if (meas_data->Gtau0_d == NULL) { return -1; }
	meas_data->G0tau_d = (double *)MKL_calloc(L*N*N, sizeof(double), MEM_DATA_ALIGN); if (meas_data->G0tau_d == NULL) { return -1; }
	meas_data->Geqlt_d = (double *)MKL_calloc(L*N*N, sizeof(double), MEM_DATA_ALIGN); if (meas_data->Geqlt_d == NULL) { return -1; }

	// density and spin correlation data for all time differences
	meas_data->nn_corr = (double *)MKL_calloc(L*Ncell*Norb*Norb, sizeof(double), MEM_DATA_ALIGN);
	meas_data->zz_corr = (double *)MKL_calloc(L*Ncell*Norb*Norb, sizeof(double), MEM_DATA_ALIGN);
	meas_data->xx_corr = (double *)MKL_calloc(L*Ncell*Norb*Norb, sizeof(double), MEM_DATA_ALIGN);

	// superconducting susceptibilities
	meas_data->sc_c_sw = (double *)MKL_calloc(L*N, sizeof(double), MEM_DATA_ALIGN);
	meas_data->sc_c_dw = (double *)MKL_calloc(L*N, sizeof(double), MEM_DATA_ALIGN);
	meas_data->sc_c_sx = (double *)MKL_calloc(L*N, sizeof(double), MEM_DATA_ALIGN);

	// construct lattice coordinate sum map
	meas_data->latt_sum_map = (int *)MKL_malloc(Ncell*Ncell * sizeof(int), MEM_DATA_ALIGN);
	ConstructLatticeCoordinateSumMap(Nx, Ny, meas_data->latt_sum_map);

	// construct lattice nearest neighbor map
	meas_data->latt_xp1_map = (int *)MKL_malloc(Ncell * sizeof(int), MEM_DATA_ALIGN);
	meas_data->latt_xm1_map = (int *)MKL_malloc(Ncell * sizeof(int), MEM_DATA_ALIGN);
	meas_data->latt_yp1_map = (int *)MKL_malloc(Ncell * sizeof(int), MEM_DATA_ALIGN);
	meas_data->latt_ym1_map = (int *)MKL_malloc(Ncell * sizeof(int), MEM_DATA_ALIGN);
	ConstructLatticeNearestNeighborMap(Nx, Ny, meas_data->latt_xp1_map, meas_data->latt_xm1_map, meas_data->latt_yp1_map, meas_data->latt_ym1_map);

	meas_data->sign = 0;

	// no samples collected so far
	meas_data->nsampl = 0;

	return 0;
}


//________________________________________________________________________________________________________________________
///
/// \brief Delete unequal time measurement data (free memory)
///
void DeleteUnequalTimeMeasurementData(measurement_data_unequal_time_t *restrict meas_data)
{
	MKL_free(meas_data->latt_ym1_map);
	MKL_free(meas_data->latt_yp1_map);
	MKL_free(meas_data->latt_xm1_map);
	MKL_free(meas_data->latt_xp1_map);

	MKL_free(meas_data->latt_sum_map);

	MKL_free(meas_data->Hd);
	MKL_free(meas_data->Hu);

	MKL_free(meas_data->sc_c_sw);
	MKL_free(meas_data->sc_c_dw);
	MKL_free(meas_data->sc_c_sx);

	MKL_free(meas_data->xx_corr);
	MKL_free(meas_data->zz_corr);
	MKL_free(meas_data->nn_corr);

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
void AccumulateUnequalTimeMeasurement(const double sign, const double *const *Bu, const double *const *Bd, measurement_data_unequal_time_t *restrict meas_data)
{
	int l;
	int i, k;
	int o, p;

	const int Norb  = meas_data->Norb;
	const int Ncell = meas_data->Ncell;
	const int N = Norb * Ncell;
	const int L = meas_data->L;

	// current spin-up unequal time Green's functions
	double *curGtau0_u = (double *)MKL_malloc(L*N*N * sizeof(double), MEM_DATA_ALIGN);
	double *curG0tau_u = (double *)MKL_malloc(L*N*N * sizeof(double), MEM_DATA_ALIGN);
	double *curGeqlt_u = (double *)MKL_malloc(L*N*N * sizeof(double), MEM_DATA_ALIGN);

	// current spin-down unequal time Green's functions
	double *curGtau0_d = (double *)MKL_malloc(L*N*N * sizeof(double), MEM_DATA_ALIGN);
	double *curG0tau_d = (double *)MKL_malloc(L*N*N * sizeof(double), MEM_DATA_ALIGN);
	double *curGeqlt_d = (double *)MKL_malloc(L*N*N * sizeof(double), MEM_DATA_ALIGN);

	// compute and accumulate unequal time spin-up and spin-down Green's functions
	#pragma omp parallel sections
	{
		#pragma omp section
		{
			ComputeUnequalTimeGreensFunction(N, L, Bu, meas_data->Hu, curGtau0_u, curG0tau_u, curGeqlt_u, NULL, NULL, NULL);
			cblas_daxpy(L*N*N, sign, curGtau0_u, 1, meas_data->Gtau0_u, 1);
			cblas_daxpy(L*N*N, sign, curG0tau_u, 1, meas_data->G0tau_u, 1);
			cblas_daxpy(L*N*N, sign, curGeqlt_u, 1, meas_data->Geqlt_u, 1);
		}
		#pragma omp section
		{
			ComputeUnequalTimeGreensFunction(N, L, Bd, meas_data->Hd, curGtau0_d, curG0tau_d, curGeqlt_d, NULL, NULL, NULL);
			cblas_daxpy(L*N*N, sign, curGtau0_d, 1, meas_data->Gtau0_d, 1);
			cblas_daxpy(L*N*N, sign, curG0tau_d, 1, meas_data->G0tau_d, 1);
			cblas_daxpy(L*N*N, sign, curGeqlt_d, 1, meas_data->Geqlt_d, 1);
		}
	}

	// accumulate density and spin correlation data

	// sign and normalization factor
	const double signfac = sign / Ncell;

	for (l = 0; l < L; l++)		// for all discrete time differences...
	{
		for (o = 0; o < Norb; o++)
		{
			for (p = 0; p < Norb; p++)
			{
				const int offset = (o + p*Norb + Norb*Norb*l) * Ncell;

				for (i = 0; i < Ncell; i++)
				{
					const int io = i + o * Ncell;
					const double Gtt_u_ii = curGeqlt_u[io + N*(io + N*l)];
					const double Gtt_d_ii = curGeqlt_d[io + N*(io + N*l)];

					for (k = 0; k < Ncell; k++)
					{
						// special case: same lattice site, orbital and time slice
						if (k == 0 && o == p && l == 0)
						{
							meas_data->nn_corr[0 + offset] += signfac*(4 - 3*(Gtt_u_ii + Gtt_d_ii) + 2*Gtt_u_ii*Gtt_d_ii);
							meas_data->zz_corr[0 + offset] += signfac*(Gtt_u_ii - 2*Gtt_u_ii*Gtt_d_ii + Gtt_d_ii);
							meas_data->xx_corr[0 + offset] += signfac*(Gtt_u_ii - 2*Gtt_u_ii*Gtt_d_ii + Gtt_d_ii);
						}
						else
						{
							const int jp = meas_data->latt_sum_map[k + Ncell*i] + p*Ncell;

							const double G00_u_jj = curGeqlt_u[jp + N*jp];
							const double G00_d_jj = curGeqlt_d[jp + N*jp];

							const double Gt0_u_ij = curGtau0_u[io + N*(L*jp + l)];
							const double Gt0_d_ij = curGtau0_d[io + N*(L*jp + l)];
							const double G0t_u_ji = curG0tau_u[jp + N*(io + N*l)];
							const double G0t_d_ji = curG0tau_d[jp + N*(io + N*l)];

							meas_data->nn_corr[k + offset] += signfac*((2 - Gtt_u_ii - Gtt_d_ii)*(2 - G00_u_jj - G00_d_jj) - (Gt0_u_ij*G0t_u_ji + Gt0_d_ij*G0t_d_ji));
							meas_data->zz_corr[k + offset] += signfac*((    Gtt_u_ii - Gtt_d_ii)*(    G00_u_jj - G00_d_jj) - (Gt0_u_ij*G0t_u_ji + Gt0_d_ij*G0t_d_ji));
							meas_data->xx_corr[k + offset] += signfac*(                                                    - (Gt0_u_ij*G0t_d_ji + Gt0_d_ij*G0t_u_ji));
						}
					}
				}
			}


			// superconducting susceptibilities (separately for each orbital type)

			// base pointer of unequal time Green's functions G(tau, 0) for current time difference tau = l and orbital o
			const double *Gt0_u_base = &curGtau0_u[(o*Ncell) + N*(L*(o*Ncell) + l)];
			const double *Gt0_d_base = &curGtau0_d[(o*Ncell) + N*(L*(o*Ncell) + l)];

			const int offset = o*Ncell + l*N;

			for (i = 0; i < Ncell; i++)
			{
				// nearest neighbors of lattice site i
				const int ipx = meas_data->latt_xp1_map[i];
				const int imx = meas_data->latt_xm1_map[i];
				const int ipy = meas_data->latt_yp1_map[i];
				const int imy = meas_data->latt_ym1_map[i];

				for (k = 0; k < Ncell; k++)
				{
					const int j = meas_data->latt_sum_map[k + Ncell*i];

					// nearest neighbors of lattice site j
					const int jpx = meas_data->latt_xp1_map[j];
					const int jmx = meas_data->latt_xm1_map[j];
					const int jpy = meas_data->latt_yp1_map[j];
					const int jmy = meas_data->latt_ym1_map[j];

					const double Gt0_u_ij = Gt0_u_base[i + N*L*j];
					const double Gt0_d_ij = Gt0_d_base[i + N*L*j];

					// s-wave
					meas_data->sc_c_sw[k + offset] += signfac*(Gt0_u_ij * Gt0_d_ij);

					// d-wave: all 16 combinations of nearest neighbors of i and j
					meas_data->sc_c_dw[k + offset] += signfac * 0.25 * Gt0_u_ij*(
						+ Gt0_d_base[ipx + N*L*jpx] + Gt0_d_base[ipx + N*L*jmx] - Gt0_d_base[ipx + N*L*jpy] - Gt0_d_base[ipx + N*L*jmy]
						+ Gt0_d_base[imx + N*L*jpx] + Gt0_d_base[imx + N*L*jmx] - Gt0_d_base[imx + N*L*jpy] - Gt0_d_base[imx + N*L*jmy]
						- Gt0_d_base[ipy + N*L*jpx] - Gt0_d_base[ipy + N*L*jmx] + Gt0_d_base[ipy + N*L*jpy] + Gt0_d_base[ipy + N*L*jmy]
						- Gt0_d_base[imy + N*L*jpx] - Gt0_d_base[imy + N*L*jmx] + Gt0_d_base[imy + N*L*jpy] + Gt0_d_base[imy + N*L*jmy]);

					// extended s-wave: similar to d-wave, but without sign flip
					meas_data->sc_c_dw[k + offset] += signfac * 0.25 * Gt0_u_ij*(
						+ Gt0_d_base[ipx + N*L*jpx] + Gt0_d_base[ipx + N*L*jmx] + Gt0_d_base[ipx + N*L*jpy] + Gt0_d_base[ipx + N*L*jmy]
						+ Gt0_d_base[imx + N*L*jpx] + Gt0_d_base[imx + N*L*jmx] + Gt0_d_base[imx + N*L*jpy] + Gt0_d_base[imx + N*L*jmy]
						+ Gt0_d_base[ipy + N*L*jpx] + Gt0_d_base[ipy + N*L*jmx] + Gt0_d_base[ipy + N*L*jpy] + Gt0_d_base[ipy + N*L*jmy]
						+ Gt0_d_base[imy + N*L*jpx] + Gt0_d_base[imy + N*L*jmx] + Gt0_d_base[imy + N*L*jpy] + Gt0_d_base[imy + N*L*jmy]);
				}
			}
		}
	}

	// clean up
	MKL_free(curGeqlt_d);
	MKL_free(curG0tau_d);
	MKL_free(curGtau0_d);
	MKL_free(curGeqlt_u);
	MKL_free(curG0tau_u);
	MKL_free(curGtau0_u);

	// add current sign
	meas_data->sign += sign;

	// increment sample counter
	meas_data->nsampl++;
}


//________________________________________________________________________________________________________________________
///
/// \brief Normalize unequal time measurement data (divide by accumulated sign)
///
void NormalizeUnequalTimeMeasurementData(measurement_data_unequal_time_t *meas_data)
{
	// total number of orbitals and cells
	const int Norb  = meas_data->Norb;
	const int Ncell = meas_data->Ncell;

	const int N = Norb * Ncell;

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

	// divide density and spin correlations by sign
	cblas_dscal(L*Ncell*Norb*Norb, nfac, meas_data->nn_corr, 1);
	cblas_dscal(L*Ncell*Norb*Norb, nfac, meas_data->zz_corr, 1);
	cblas_dscal(L*Ncell*Norb*Norb, nfac, meas_data->xx_corr, 1);

	// divide superconducting susceptibilities by sign
	cblas_dscal(L*N, nfac, meas_data->sc_c_sw, 1);
	cblas_dscal(L*N, nfac, meas_data->sc_c_dw, 1);
	cblas_dscal(L*N, nfac, meas_data->sc_c_sx, 1);

	// calculate average sign
	meas_data->sign /= meas_data->nsampl;
}


//________________________________________________________________________________________________________________________
///
/// \brief Load unequal time measurement data structure from disk
///
void LoadUnequalTimeMeasurementData(const char *fnbase, const measurement_data_unequal_time_t *meas_data)
{
	const int Norb  = meas_data->Norb;
	const int Ncell = meas_data->Ncell;
	const int N = Norb * Ncell;
	const int L = meas_data->L;

	char path[1024];
	sprintf(path, "%s_uneqlt_sign.dat",    fnbase); ReadData(path, (void *)&meas_data->sign,   sizeof(double), 1);
	sprintf(path, "%s_uneqlt_nsampl.dat",  fnbase); ReadData(path, (void *)&meas_data->nsampl, sizeof(int), 1);

	sprintf(path, "%s_uneqlt_Gtau0_u.dat", fnbase); ReadData(path, meas_data->Gtau0_u, sizeof(double), L*N*N);
	sprintf(path, "%s_uneqlt_G0tau_u.dat", fnbase); ReadData(path, meas_data->G0tau_u, sizeof(double), L*N*N);
	sprintf(path, "%s_uneqlt_Geqlt_u.dat", fnbase); ReadData(path, meas_data->Geqlt_u, sizeof(double), L*N*N);
	sprintf(path, "%s_uneqlt_Gtau0_d.dat", fnbase); ReadData(path, meas_data->Gtau0_d, sizeof(double), L*N*N);
	sprintf(path, "%s_uneqlt_G0tau_d.dat", fnbase); ReadData(path, meas_data->G0tau_d, sizeof(double), L*N*N);
	sprintf(path, "%s_uneqlt_Geqlt_d.dat", fnbase); ReadData(path, meas_data->Geqlt_d, sizeof(double), L*N*N);

	sprintf(path, "%s_uneqlt_nn_corr.dat", fnbase); ReadData(path, meas_data->nn_corr, sizeof(double), L*Ncell*Norb*Norb);
	sprintf(path, "%s_uneqlt_zz_corr.dat", fnbase); ReadData(path, meas_data->zz_corr, sizeof(double), L*Ncell*Norb*Norb);
	sprintf(path, "%s_uneqlt_xx_corr.dat", fnbase); ReadData(path, meas_data->xx_corr, sizeof(double), L*Ncell*Norb*Norb);

	sprintf(path, "%s_uneqlt_sc_c_sw.dat", fnbase); ReadData(path, meas_data->sc_c_sw, sizeof(double), L*N);
	sprintf(path, "%s_uneqlt_sc_c_dw.dat", fnbase); ReadData(path, meas_data->sc_c_dw, sizeof(double), L*N);
	sprintf(path, "%s_uneqlt_sc_c_sx.dat", fnbase); ReadData(path, meas_data->sc_c_sx, sizeof(double), L*N);
}


//________________________________________________________________________________________________________________________
///
/// \brief Save unequal time measurement data structure to disk
///
void SaveUnequalTimeMeasurementData(const char *fnbase, const measurement_data_unequal_time_t *meas_data)
{
	const int Norb  = meas_data->Norb;
	const int Ncell = meas_data->Ncell;
	const int N = Norb * Ncell;
	const int L = meas_data->L;

	char path[1024];
	sprintf(path, "%s_uneqlt_sign.dat",    fnbase); WriteData(path, &meas_data->sign,   sizeof(double), 1, false);
	sprintf(path, "%s_uneqlt_nsampl.dat",  fnbase); WriteData(path, &meas_data->nsampl, sizeof(int),    1, false);

	sprintf(path, "%s_uneqlt_Gtau0_u.dat", fnbase); WriteData(path, meas_data->Gtau0_u, sizeof(double), L*N*N, false);
	sprintf(path, "%s_uneqlt_G0tau_u.dat", fnbase); WriteData(path, meas_data->G0tau_u, sizeof(double), L*N*N, false);
	sprintf(path, "%s_uneqlt_Geqlt_u.dat", fnbase); WriteData(path, meas_data->Geqlt_u, sizeof(double), L*N*N, false);
	sprintf(path, "%s_uneqlt_Gtau0_d.dat", fnbase); WriteData(path, meas_data->Gtau0_d, sizeof(double), L*N*N, false);
	sprintf(path, "%s_uneqlt_G0tau_d.dat", fnbase); WriteData(path, meas_data->G0tau_d, sizeof(double), L*N*N, false);
	sprintf(path, "%s_uneqlt_Geqlt_d.dat", fnbase); WriteData(path, meas_data->Geqlt_d, sizeof(double), L*N*N, false);

	sprintf(path, "%s_uneqlt_nn_corr.dat", fnbase); WriteData(path, meas_data->nn_corr, sizeof(double), L*Ncell*Norb*Norb, false);
	sprintf(path, "%s_uneqlt_zz_corr.dat", fnbase); WriteData(path, meas_data->zz_corr, sizeof(double), L*Ncell*Norb*Norb, false);
	sprintf(path, "%s_uneqlt_xx_corr.dat", fnbase); WriteData(path, meas_data->xx_corr, sizeof(double), L*Ncell*Norb*Norb, false);

	sprintf(path, "%s_uneqlt_sc_c_sw.dat", fnbase); WriteData(path, meas_data->sc_c_sw, sizeof(double), L*N, false);
	sprintf(path, "%s_uneqlt_sc_c_dw.dat", fnbase); WriteData(path, meas_data->sc_c_dw, sizeof(double), L*N, false);
	sprintf(path, "%s_uneqlt_sc_c_sx.dat", fnbase); WriteData(path, meas_data->sc_c_sx, sizeof(double), L*N, false);
}
