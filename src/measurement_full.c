// measurement.c assumes translational symmetry, which is appropriate when the
// lattice has periodic boundary conditions. If you're using open or
// cylindrical boundary conditions or doing something else that breaks
// translational symmetry, replace measurement.c with this file in order to
// record measurement data for every site.
//
// The d-wave superconducting susceptibility still relies on a rectangular lattice geometry.

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
/// \brief Construct nearest neighbor maps for a rectangular lattice
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
void AllocateMeasurementData(const int Norb, const int Nx, const int Ny, const int pbc_shift, measurement_data_t *restrict meas_data)
{
	// unused parameter
	(void)pbc_shift;

	const int Ncell = Nx * Ny;
	const int N     = Ncell * Norb;

	meas_data->Norb  = Norb;
	meas_data->Ncell = Ncell;

	// no samples collected so far
	meas_data->nsampl = 0;

	meas_data->sign = 0;

	meas_data->density_u = (double *)MKL_calloc(N, sizeof(double), MEM_DATA_ALIGN);
	meas_data->density_d = (double *)MKL_calloc(N, sizeof(double), MEM_DATA_ALIGN);
	meas_data->doubleocc = (double *)MKL_calloc(N, sizeof(double), MEM_DATA_ALIGN);

	// green's functions
	meas_data->grfun_u = (double *)MKL_calloc(N*N, sizeof(double), MEM_DATA_ALIGN);
	meas_data->grfun_d = (double *)MKL_calloc(N*N, sizeof(double), MEM_DATA_ALIGN);

	// density correlation data
	meas_data->uu_corr = (double *)MKL_calloc(N*N, sizeof(double), MEM_DATA_ALIGN);
	meas_data->dd_corr = (double *)MKL_calloc(N*N, sizeof(double), MEM_DATA_ALIGN);
	meas_data->ud_corr = (double *)MKL_calloc(N*N, sizeof(double), MEM_DATA_ALIGN);
	meas_data->ff_corr = (double *)MKL_calloc(N*N, sizeof(double), MEM_DATA_ALIGN);

	// spin correlation data
	meas_data->zz_corr = (double *)MKL_calloc(N*N, sizeof(double), MEM_DATA_ALIGN);
	meas_data->xx_corr = (double *)MKL_calloc(N*N, sizeof(double), MEM_DATA_ALIGN);

	// superconducting susceptibilities
	meas_data->sc_c_sw = (double *)MKL_calloc(Ncell*Ncell*Norb, sizeof(double), MEM_DATA_ALIGN);
	meas_data->sc_c_dw = (double *)MKL_calloc(Ncell*Ncell*Norb, sizeof(double), MEM_DATA_ALIGN);
	meas_data->sc_c_sx = (double *)MKL_calloc(Ncell*Ncell*Norb, sizeof(double), MEM_DATA_ALIGN);

	// construct lattice nearest neighbor map
	meas_data->latt_xp1_map = (int *)MKL_malloc(Ncell * sizeof(int), MEM_DATA_ALIGN);
	meas_data->latt_xm1_map = (int *)MKL_malloc(Ncell * sizeof(int), MEM_DATA_ALIGN);
	meas_data->latt_yp1_map = (int *)MKL_malloc(Ncell * sizeof(int), MEM_DATA_ALIGN);
	meas_data->latt_ym1_map = (int *)MKL_malloc(Ncell * sizeof(int), MEM_DATA_ALIGN);
	ConstructLatticeNearestNeighborMap(Nx, Ny, meas_data->latt_xp1_map, meas_data->latt_xm1_map, meas_data->latt_yp1_map, meas_data->latt_ym1_map);
}


//________________________________________________________________________________________________________________________
///
/// \brief Delete measurement data (free memory)
///
void DeleteMeasurementData(measurement_data_t *restrict meas_data)
{
	MKL_free(meas_data->latt_ym1_map);
	MKL_free(meas_data->latt_yp1_map);
	MKL_free(meas_data->latt_xm1_map);
	MKL_free(meas_data->latt_xp1_map);

	MKL_free(meas_data->sc_c_sw);
	MKL_free(meas_data->sc_c_dw);
	MKL_free(meas_data->sc_c_sx);

	MKL_free(meas_data->xx_corr);
	MKL_free(meas_data->zz_corr);
	MKL_free(meas_data->ff_corr);
	MKL_free(meas_data->ud_corr);
	MKL_free(meas_data->dd_corr);
	MKL_free(meas_data->uu_corr);

	MKL_free(meas_data->grfun_d);
	MKL_free(meas_data->grfun_u);

	MKL_free(meas_data->doubleocc);
	MKL_free(meas_data->density_d);
	MKL_free(meas_data->density_u);
}


//________________________________________________________________________________________________________________________
///
/// \brief Accumulate equal time "measurement" data
///
void AccumulateMeasurement(const greens_func_t *restrict Gu, const greens_func_t *restrict Gd, measurement_data_t *restrict meas_data)
{
	int i, j;

	// lattice dimensions
	const int Norb  = meas_data->Norb;
	const int Ncell = meas_data->Ncell;
	const int N     = Ncell * Norb;

	// product of the determinant signs of the Green's function matrices
	const double sign = (double)(Gu->sgndet * Gd->sgndet);
	assert(sign != 0);

	for (i = 0; i < N; i++)
	{
		meas_data->density_u[i] += sign*(1 - Gu->mat[i + N*i]);
		meas_data->density_d[i] += sign*(1 - Gd->mat[i + N*i]);
		meas_data->doubleocc[i] += sign*(1 - Gu->mat[i + N*i])*(1 - Gd->mat[i + N*i]);
	}

	// green's functions
	cblas_daxpy(N*N, sign, Gu->mat, 1, meas_data->grfun_u, 1);
	cblas_daxpy(N*N, sign, Gd->mat, 1, meas_data->grfun_d, 1);

	// density and spin correlations

	for (i = 0; i < N; i++)
	{
		const double Gu_ii = Gu->mat[i + N*i];
		const double Gd_ii = Gd->mat[i + N*i];

		for (j = 0; j < N; j++)
		{
			const double Gu_jj = Gu->mat[j + N*j];
			const double Gd_jj = Gd->mat[j + N*j];
			const double Gu_ij = Gu->mat[i + N*j];
			const double Gd_ij = Gd->mat[i + N*j];
			const double Gu_ji = Gu->mat[j + N*i];
			const double Gd_ji = Gd->mat[j + N*i];

			if (i == j)		// special case: same lattice site and orbital
			{
				meas_data->uu_corr[i + N*j] += sign*(1 - Gu_ii);
				meas_data->dd_corr[i + N*j] += sign*(1 - Gd_ii);
				meas_data->ff_corr[i + N*j] += sign*(Gu_ii*(1 - Gu_ii) + Gd_ii*(1 - Gd_ii));
				meas_data->zz_corr[i + N*j] += sign*(Gu_ii - 2*Gu_ii*Gd_ii + Gd_ii);
				meas_data->xx_corr[i + N*j] += sign*(Gu_ii - 2*Gu_ii*Gd_ii + Gd_ii);
			}
			else
			{
				meas_data->uu_corr[i + N*j] += sign*((1 - Gu_ii)*(1 - Gu_jj) - Gu_ij*Gu_ji);
				meas_data->dd_corr[i + N*j] += sign*((1 - Gd_ii)*(1 - Gd_jj) - Gd_ij*Gd_ji);
				meas_data->ff_corr[i + N*j] += sign*(                                - (Gu_ij*Gu_ji + Gd_ij*Gd_ji));
				meas_data->zz_corr[i + N*j] += sign*((Gu_ii - Gd_ii)*(Gu_jj - Gd_jj) - (Gu_ij*Gu_ji + Gd_ij*Gd_ji));
				meas_data->xx_corr[i + N*j] += sign*(                                - (Gu_ij*Gd_ji + Gd_ij*Gu_ji));
			}

			// independent of special case
			meas_data->ud_corr[i + N*j] += sign*((1 - Gu_ii)*(1 - Gd_jj) + (1 - Gd_ii)*(1 - Gu_jj));
		}
	}

	// superconducting susceptibilities (separately for each orbital type)

	int o;
	for (o = 0; o < Norb; o++)
	{
		// base pointer of Green's functions for current orbital type
		const double *Gu_orb = &Gu->mat[(Ncell*o) + N*(Ncell*o)];
		const double *Gd_orb = &Gd->mat[(Ncell*o) + N*(Ncell*o)];

		for (i = 0; i < Ncell; i++)
		{
			// nearest neighbors of lattice site i
			const int ipx = meas_data->latt_xp1_map[i];
			const int imx = meas_data->latt_xm1_map[i];
			const int ipy = meas_data->latt_yp1_map[i];
			const int imy = meas_data->latt_ym1_map[i];

			for (j = 0; j < Ncell; j++)
			{
				// nearest neighbors of lattice site j
				const int jpx = meas_data->latt_xp1_map[j];
				const int jmx = meas_data->latt_xm1_map[j];
				const int jpy = meas_data->latt_yp1_map[j];
				const int jmy = meas_data->latt_ym1_map[j];

				const double Gu_ij = Gu_orb[i + N*j];
				const double Gd_ij = Gd_orb[i + N*j];

				// s-wave
				meas_data->sc_c_sw[i + Ncell*j + Ncell*Ncell*o] += sign*(Gu_ij * Gd_ij);

				// d-wave: all 16 combinations of nearest neighbors of i and j
				meas_data->sc_c_dw[i + Ncell*j + Ncell*Ncell*o] += sign * 0.25 * Gu_ij*(
					+ Gd_orb[ipx + N*jpx] + Gd_orb[ipx + N*jmx] - Gd_orb[ipx + N*jpy] - Gd_orb[ipx + N*jmy]
					+ Gd_orb[imx + N*jpx] + Gd_orb[imx + N*jmx] - Gd_orb[imx + N*jpy] - Gd_orb[imx + N*jmy]
					- Gd_orb[ipy + N*jpx] - Gd_orb[ipy + N*jmx] + Gd_orb[ipy + N*jpy] + Gd_orb[ipy + N*jmy]
					- Gd_orb[imy + N*jpx] - Gd_orb[imy + N*jmx] + Gd_orb[imy + N*jpy] + Gd_orb[imy + N*jmy]);

				// extended s-wave: similar to d-wave, but without sign flip
				meas_data->sc_c_sx[i + Ncell*j + Ncell*Ncell*o] += sign * 0.25 * Gu_ij*(
					+ Gd_orb[ipx + N*jpx] + Gd_orb[ipx + N*jmx] + Gd_orb[ipx + N*jpy] + Gd_orb[ipx + N*jmy]
					+ Gd_orb[imx + N*jpx] + Gd_orb[imx + N*jmx] + Gd_orb[imx + N*jpy] + Gd_orb[imx + N*jmy]
					+ Gd_orb[ipy + N*jpx] + Gd_orb[ipy + N*jmx] + Gd_orb[ipy + N*jpy] + Gd_orb[ipy + N*jmy]
					+ Gd_orb[imy + N*jpx] + Gd_orb[imy + N*jmx] + Gd_orb[imy + N*jpy] + Gd_orb[imy + N*jmy]);
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
/// \brief Normalize measurement data (divide by number of samples)
///
void NormalizeMeasurementData(measurement_data_t *meas_data)
{
	// total number of orbitals and cells
	const int Norb  = meas_data->Norb;
	const int Ncell = meas_data->Ncell;
	const int N     = Ncell * Norb;

	// normalization factor
	const double nfac = 1.0 / meas_data->nsampl;

	// divide by nsampl to get <s A>
	cblas_dscal(N, nfac, meas_data->density_u, 1);
	cblas_dscal(N, nfac, meas_data->density_d, 1);
	cblas_dscal(N, nfac, meas_data->doubleocc, 1);
	cblas_dscal(N*N, nfac, meas_data->grfun_u, 1);
	cblas_dscal(N*N, nfac, meas_data->grfun_d, 1);
	int m = N*N;
	cblas_dscal(m, nfac, meas_data->uu_corr, 1);
	cblas_dscal(m, nfac, meas_data->dd_corr, 1);
	cblas_dscal(m, nfac, meas_data->ud_corr, 1);
	cblas_dscal(m, nfac, meas_data->ff_corr, 1);
	cblas_dscal(m, nfac, meas_data->zz_corr, 1);
	cblas_dscal(m, nfac, meas_data->xx_corr, 1);
	// superconducting susceptibilities
	m = Ncell*Ncell*Norb;
	cblas_dscal(m, nfac, meas_data->sc_c_sw, 1);
	cblas_dscal(m, nfac, meas_data->sc_c_dw, 1);
	cblas_dscal(m, nfac, meas_data->sc_c_sx, 1);

	// calculate average sign
	meas_data->sign *= nfac;
}


//________________________________________________________________________________________________________________________
///
/// \brief Print a summary of basic measurements
///
void PrintMeasurementDataSummary(const measurement_data_t *meas_data)
{
	duprintf("_______________________________________________________________________________\n");
	duprintf("Summary of simulation results\n\n");
	duprintf("                                  average sign: %g\n", meas_data->sign);
	double total_density = 0;
	int i;
	for (i = 0; i < meas_data->Norb * meas_data->Ncell; i++)
	{
		total_density += meas_data->density_u[i] + meas_data->density_d[i];
	}
	total_density /= meas_data->Ncell;
	duprintf("         average total density (per unit cell): %g\n", total_density);
}


//________________________________________________________________________________________________________________________
///
/// \brief Load equal time measurement data structure from disk
///
void LoadMeasurementData(const char *fnbase, measurement_data_t *meas_data)
{
	const int Norb  = meas_data->Norb;
	const int Ncell = meas_data->Ncell;
	const int N     = Ncell * Norb;

	char path[1024];
	sprintf(path, "%s_eqlt_sign.dat",      fnbase); ReadData(path, (void *)&meas_data->sign,   sizeof(double), 1);
	sprintf(path, "%s_eqlt_nsampl.dat",    fnbase); ReadData(path, (void *)&meas_data->nsampl, sizeof(int),    1);
	sprintf(path, "%s_eqlt_density_u.dat", fnbase); ReadData(path, meas_data->density_u, sizeof(double), N);
	sprintf(path, "%s_eqlt_density_d.dat", fnbase); ReadData(path, meas_data->density_d, sizeof(double), N);
	sprintf(path, "%s_eqlt_doubleocc.dat", fnbase); ReadData(path, meas_data->doubleocc, sizeof(double), N);
	sprintf(path, "%s_eqlt_grfun_u.dat",   fnbase); ReadData(path, meas_data->grfun_u,   sizeof(double), N*N);
	sprintf(path, "%s_eqlt_grfun_d.dat",   fnbase); ReadData(path, meas_data->grfun_d,   sizeof(double), N*N);
	sprintf(path, "%s_eqlt_uu_corr.dat",   fnbase); ReadData(path, meas_data->uu_corr,   sizeof(double), N*N);
	sprintf(path, "%s_eqlt_dd_corr.dat",   fnbase); ReadData(path, meas_data->dd_corr,   sizeof(double), N*N);
	sprintf(path, "%s_eqlt_ud_corr.dat",   fnbase); ReadData(path, meas_data->ud_corr,   sizeof(double), N*N);
	sprintf(path, "%s_eqlt_ff_corr.dat",   fnbase); ReadData(path, meas_data->ff_corr,   sizeof(double), N*N);
	sprintf(path, "%s_eqlt_zz_corr.dat",   fnbase); ReadData(path, meas_data->zz_corr,   sizeof(double), N*N);
	sprintf(path, "%s_eqlt_xx_corr.dat",   fnbase); ReadData(path, meas_data->xx_corr,   sizeof(double), N*N);
	sprintf(path, "%s_eqlt_sc_c_sw.dat",   fnbase); ReadData(path, meas_data->sc_c_sw,   sizeof(double), Ncell*Ncell*Norb);
	sprintf(path, "%s_eqlt_sc_c_dw.dat",   fnbase); ReadData(path, meas_data->sc_c_dw,   sizeof(double), Ncell*Ncell*Norb);
	sprintf(path, "%s_eqlt_sc_c_sx.dat",   fnbase); ReadData(path, meas_data->sc_c_sx,   sizeof(double), Ncell*Ncell*Norb);
}


//________________________________________________________________________________________________________________________
///
/// \brief Save equal time measurement data structure to disk
///
void SaveMeasurementData(const char *fnbase, const measurement_data_t *meas_data)
{
	const int Norb  = meas_data->Norb;
	const int Ncell = meas_data->Ncell;
	const int N     = Ncell * Norb;

	char path[1024];
	sprintf(path, "%s_eqlt_sign.dat",      fnbase); WriteData(path, &meas_data->sign,     sizeof(double), 1, false);
	sprintf(path, "%s_eqlt_nsampl.dat",    fnbase); WriteData(path, &meas_data->nsampl,   sizeof(int),    1, false);
	sprintf(path, "%s_eqlt_density_u.dat", fnbase); WriteData(path, meas_data->density_u, sizeof(double), N, false);
	sprintf(path, "%s_eqlt_density_d.dat", fnbase); WriteData(path, meas_data->density_d, sizeof(double), N, false);
	sprintf(path, "%s_eqlt_doubleocc.dat", fnbase); WriteData(path, meas_data->doubleocc, sizeof(double), N, false);
	sprintf(path, "%s_eqlt_grfun_u.dat",   fnbase); WriteData(path, meas_data->grfun_u,   sizeof(double), N*N, false);
	sprintf(path, "%s_eqlt_grfun_d.dat",   fnbase); WriteData(path, meas_data->grfun_d,   sizeof(double), N*N, false);
	sprintf(path, "%s_eqlt_uu_corr.dat",   fnbase); WriteData(path, meas_data->uu_corr,   sizeof(double), N*N, false);
	sprintf(path, "%s_eqlt_dd_corr.dat",   fnbase); WriteData(path, meas_data->dd_corr,   sizeof(double), N*N, false);
	sprintf(path, "%s_eqlt_ud_corr.dat",   fnbase); WriteData(path, meas_data->ud_corr,   sizeof(double), N*N, false);
	sprintf(path, "%s_eqlt_ff_corr.dat",   fnbase); WriteData(path, meas_data->ff_corr,   sizeof(double), N*N, false);
	sprintf(path, "%s_eqlt_zz_corr.dat",   fnbase); WriteData(path, meas_data->zz_corr,   sizeof(double), N*N, false);
	sprintf(path, "%s_eqlt_xx_corr.dat",   fnbase); WriteData(path, meas_data->xx_corr,   sizeof(double), N*N, false);
	sprintf(path, "%s_eqlt_sc_c_sw.dat",   fnbase); WriteData(path, meas_data->sc_c_sw,   sizeof(double), Ncell*Ncell*Norb, false);
	sprintf(path, "%s_eqlt_sc_c_dw.dat",   fnbase); WriteData(path, meas_data->sc_c_dw,   sizeof(double), Ncell*Ncell*Norb, false);
	sprintf(path, "%s_eqlt_sc_c_sx.dat",   fnbase); WriteData(path, meas_data->sc_c_sx,   sizeof(double), Ncell*Ncell*Norb, false);
}


//________________________________________________________________________________________________________________________
///
/// \brief Allocate and initialize unequal time measurement data structure
///
int AllocateUnequalTimeMeasurementData(const int Norb, const int Nx, const int Ny, const int pbc_shift, const int L, const int numBprod, measurement_data_unequal_time_t *restrict meas_data)
{
	// unused parameter
	(void)pbc_shift;

	// lattice dimensions
	const int Ncell = Nx * Ny;
	const int N     = Ncell * Norb;

	meas_data->Norb  = Norb;
	meas_data->Ncell = Ncell;

	// total number of time steps
	meas_data->L = L;

	// allocate temporary H matrices
	meas_data->Hu = (double *)MKL_malloc(N*N*numBprod*numBprod * sizeof(double), MEM_DATA_ALIGN); if (meas_data->Hu == NULL) { return -1; }
	meas_data->Hd = (double *)MKL_malloc(N*N*numBprod*numBprod * sizeof(double), MEM_DATA_ALIGN); if (meas_data->Hd == NULL) { return -1; }

	// allocate and initialize Green's functions with zeros
	meas_data->Gtau0_u = (double *)MKL_calloc(N*L*N, sizeof(double), MEM_DATA_ALIGN); if (meas_data->Gtau0_u == NULL) { return -1; }
	meas_data->G0tau_u = (double *)MKL_calloc(N*N*L, sizeof(double), MEM_DATA_ALIGN); if (meas_data->G0tau_u == NULL) { return -1; }
	meas_data->Geqlt_u = (double *)MKL_calloc(N*N*L, sizeof(double), MEM_DATA_ALIGN); if (meas_data->Geqlt_u == NULL) { return -1; }
	meas_data->Gtau0_d = (double *)MKL_calloc(N*L*N, sizeof(double), MEM_DATA_ALIGN); if (meas_data->Gtau0_d == NULL) { return -1; }
	meas_data->G0tau_d = (double *)MKL_calloc(N*N*L, sizeof(double), MEM_DATA_ALIGN); if (meas_data->G0tau_d == NULL) { return -1; }
	meas_data->Geqlt_d = (double *)MKL_calloc(N*N*L, sizeof(double), MEM_DATA_ALIGN); if (meas_data->Geqlt_d == NULL) { return -1; }

	// density and spin correlation data for all time differences
	meas_data->nn_corr = (double *)MKL_calloc(N*N*L, sizeof(double), MEM_DATA_ALIGN);
	meas_data->zz_corr = (double *)MKL_calloc(N*N*L, sizeof(double), MEM_DATA_ALIGN);
	meas_data->xx_corr = (double *)MKL_calloc(N*N*L, sizeof(double), MEM_DATA_ALIGN);

	// superconducting susceptibilities
	meas_data->sc_c_sw = (double *)MKL_calloc(Ncell*Ncell*Norb*L, sizeof(double), MEM_DATA_ALIGN);
	meas_data->sc_c_dw = (double *)MKL_calloc(Ncell*Ncell*Norb*L, sizeof(double), MEM_DATA_ALIGN);
	meas_data->sc_c_sx = (double *)MKL_calloc(Ncell*Ncell*Norb*L, sizeof(double), MEM_DATA_ALIGN);

	// Raman
	meas_data->ram_b1g = (double *)MKL_calloc(Ncell*Ncell*Norb*L, sizeof(double), MEM_DATA_ALIGN);
	meas_data->ram_b2g = (double *)MKL_calloc(Ncell*Ncell*Norb*L, sizeof(double), MEM_DATA_ALIGN);

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

	MKL_free(meas_data->Hd);
	MKL_free(meas_data->Hu);

	MKL_free(meas_data->ram_b2g);
	MKL_free(meas_data->ram_b1g);

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
void AccumulateUnequalTimeMeasurement(const double sign, const time_step_matrices_t *restrict tsm_u, const time_step_matrices_t *restrict tsm_d, measurement_data_unequal_time_t *restrict meas_data)
{
	int l;
	int i, j;

	const int Norb  = meas_data->Norb;
	const int Ncell = meas_data->Ncell;
	const int N = Ncell * Norb;
	const int L = meas_data->L;

	// current spin-up unequal time Green's functions
	double *curGtau0_u = (double *)MKL_malloc(N*L*N * sizeof(double), MEM_DATA_ALIGN);
	double *curG0tau_u = (double *)MKL_malloc(N*N*L * sizeof(double), MEM_DATA_ALIGN);
	double *curGeqlt_u = (double *)MKL_malloc(N*N*L * sizeof(double), MEM_DATA_ALIGN);

	// current spin-down unequal time Green's functions
	double *curGtau0_d = (double *)MKL_malloc(N*L*N * sizeof(double), MEM_DATA_ALIGN);
	double *curG0tau_d = (double *)MKL_malloc(N*N*L * sizeof(double), MEM_DATA_ALIGN);
	double *curGeqlt_d = (double *)MKL_malloc(N*N*L * sizeof(double), MEM_DATA_ALIGN);

	// compute and accumulate unequal time spin-up and spin-down Green's functions
	#pragma omp parallel sections
	{
		#pragma omp section
		{
			ComputeUnequalTimeGreensFunction(N, L, tsm_u, meas_data->Hu, curGtau0_u, curG0tau_u, curGeqlt_u);
			cblas_daxpy(N*L*N, sign, curGtau0_u, 1, meas_data->Gtau0_u, 1);
			cblas_daxpy(N*N*L, sign, curG0tau_u, 1, meas_data->G0tau_u, 1);
			cblas_daxpy(N*N*L, sign, curGeqlt_u, 1, meas_data->Geqlt_u, 1);
		}
		#pragma omp section
		{
			ComputeUnequalTimeGreensFunction(N, L, tsm_d, meas_data->Hd, curGtau0_d, curG0tau_d, curGeqlt_d);
			cblas_daxpy(N*L*N, sign, curGtau0_d, 1, meas_data->Gtau0_d, 1);
			cblas_daxpy(N*N*L, sign, curG0tau_d, 1, meas_data->G0tau_d, 1);
			cblas_daxpy(N*N*L, sign, curGeqlt_d, 1, meas_data->Geqlt_d, 1);
		}
	}

	for (l = 0; l < L; l++)		// for all discrete time differences...
	{
		// accumulate density and spin correlation data

		for (i = 0; i < N; i++)
		{
			const double Gtt_u_ii = curGeqlt_u[i + N*(i + N*l)];
			const double Gtt_d_ii = curGeqlt_d[i + N*(i + N*l)];

			for (j = 0; j < N; j++)
			{
				// special case: same lattice site, orbital and time slice
				if (i == j && l == 0)
				{
					meas_data->nn_corr[i + N*j] += sign*(4 - 3*(Gtt_u_ii + Gtt_d_ii) + 2*Gtt_u_ii*Gtt_d_ii);
					meas_data->zz_corr[i + N*j] += sign*(Gtt_u_ii - 2*Gtt_u_ii*Gtt_d_ii + Gtt_d_ii);
					meas_data->xx_corr[i + N*j] += sign*(Gtt_u_ii - 2*Gtt_u_ii*Gtt_d_ii + Gtt_d_ii);
				}
				else
				{
					const double G00_u_jj = curGeqlt_u[j + N*j];
					const double G00_d_jj = curGeqlt_d[j + N*j];

					const double Gt0_u_ij = curGtau0_u[i + N*(l + L*j)];
					const double Gt0_d_ij = curGtau0_d[i + N*(l + L*j)];
					const double G0t_u_ji = curG0tau_u[j + N*(i + N*l)];
					const double G0t_d_ji = curG0tau_d[j + N*(i + N*l)];

					meas_data->nn_corr[i + N*j + N*N*l] += sign*((2 - Gtt_u_ii - Gtt_d_ii)*(2 - G00_u_jj - G00_d_jj) - (Gt0_u_ij*G0t_u_ji + Gt0_d_ij*G0t_d_ji));
					meas_data->zz_corr[i + N*j + N*N*l] += sign*((    Gtt_u_ii - Gtt_d_ii)*(    G00_u_jj - G00_d_jj) - (Gt0_u_ij*G0t_u_ji + Gt0_d_ij*G0t_d_ji));
					meas_data->xx_corr[i + N*j + N*N*l] += sign*(                                                    - (Gt0_u_ij*G0t_d_ji + Gt0_d_ij*G0t_u_ji));
				}
			}
		}

		// superconducting susceptibilities and Raman correlation functions (separately for each orbital type)

		int o;
		for (o = 0; o < Norb; o++)
		{
			// base pointer of unequal time Green's functions G(tau, 0) for current time difference tau = l and orbital o
			const double *Gt0_u_base = &curGtau0_u[(Ncell*o) + N*(l + L*(Ncell*o))];
			const double *Gt0_d_base = &curGtau0_d[(Ncell*o) + N*(l + L*(Ncell*o))];

			// base pointer of unequal time Green's functions G(0, tau) for current time difference tau = l and orbital o
			const double *G0t_u_base = &curG0tau_u[(Ncell*o) + N*((Ncell*o) + N*l)];
			const double *G0t_d_base = &curG0tau_d[(Ncell*o) + N*((Ncell*o) + N*l)];

			// base pointer of equal time Green's functions G(0, 0) for current orbital o
			const double *G00_u_base = &curGeqlt_u[(Ncell*o) + N*(Ncell*o)];
			const double *G00_d_base = &curGeqlt_d[(Ncell*o) + N*(Ncell*o)];

			// base pointer of equal time Green's functions G(tau, tau) for current time tau = l and orbital o
			const double *Gtt_u_base = &curGeqlt_u[(Ncell*o) + N*((Ncell*o) + N*l)];
			const double *Gtt_d_base = &curGeqlt_d[(Ncell*o) + N*((Ncell*o) + N*l)];

			const int offset = Ncell*Ncell * (o + Norb*l);

			for (i = 0; i < Ncell; i++)
			{
				// nearest neighbors of lattice site i
				const int ipx = meas_data->latt_xp1_map[i];
				const int imx = meas_data->latt_xm1_map[i];
				const int ipy = meas_data->latt_yp1_map[i];
				const int imy = meas_data->latt_ym1_map[i];

				// next-nearest neighbors of lattice site i
				const int ipxpy = meas_data->latt_yp1_map[ipx];
				const int ipxmy = meas_data->latt_ym1_map[ipx];
				const int imxpy = meas_data->latt_yp1_map[imx];
				const int imxmy = meas_data->latt_ym1_map[imx];

				for (j = 0; j < Ncell; j++)
				{
					// nearest neighbors of lattice site j
					const int jpx = meas_data->latt_xp1_map[j];
					const int jmx = meas_data->latt_xm1_map[j];
					const int jpy = meas_data->latt_yp1_map[j];
					const int jmy = meas_data->latt_ym1_map[j];

					// next-nearest neighbors of lattice site j
					const int jpxpy = meas_data->latt_yp1_map[jpx];
					const int jpxmy = meas_data->latt_ym1_map[jpx];
					const int jmxpy = meas_data->latt_yp1_map[jmx];
					const int jmxmy = meas_data->latt_ym1_map[jmx];

					const double Gt0_u_ij = Gt0_u_base[i + N*L*j];
					const double Gt0_d_ij = Gt0_d_base[i + N*L*j];

					// s-wave
					meas_data->sc_c_sw[i + Ncell*j + offset] += sign*(Gt0_u_ij * Gt0_d_ij);

					// d-wave: all 16 combinations of nearest neighbors of i and j
					meas_data->sc_c_dw[i + Ncell*j + offset] += sign * 0.25 * Gt0_u_ij*(
						+ Gt0_d_base[ipx + N*L*jpx] + Gt0_d_base[ipx + N*L*jmx] - Gt0_d_base[ipx + N*L*jpy] - Gt0_d_base[ipx + N*L*jmy]
						+ Gt0_d_base[imx + N*L*jpx] + Gt0_d_base[imx + N*L*jmx] - Gt0_d_base[imx + N*L*jpy] - Gt0_d_base[imx + N*L*jmy]
						- Gt0_d_base[ipy + N*L*jpx] - Gt0_d_base[ipy + N*L*jmx] + Gt0_d_base[ipy + N*L*jpy] + Gt0_d_base[ipy + N*L*jmy]
						- Gt0_d_base[imy + N*L*jpx] - Gt0_d_base[imy + N*L*jmx] + Gt0_d_base[imy + N*L*jpy] + Gt0_d_base[imy + N*L*jmy]);

					// extended s-wave: similar to d-wave, but without sign flip
					meas_data->sc_c_sx[i + Ncell*j + offset] += sign * 0.25 * Gt0_u_ij*(
						+ Gt0_d_base[ipx + N*L*jpx] + Gt0_d_base[ipx + N*L*jmx] + Gt0_d_base[ipx + N*L*jpy] + Gt0_d_base[ipx + N*L*jmy]
						+ Gt0_d_base[imx + N*L*jpx] + Gt0_d_base[imx + N*L*jmx] + Gt0_d_base[imx + N*L*jpy] + Gt0_d_base[imx + N*L*jmy]
						+ Gt0_d_base[ipy + N*L*jpx] + Gt0_d_base[ipy + N*L*jmx] + Gt0_d_base[ipy + N*L*jpy] + Gt0_d_base[ipy + N*L*jmy]
						+ Gt0_d_base[imy + N*L*jpx] + Gt0_d_base[imy + N*L*jmx] + Gt0_d_base[imy + N*L*jpy] + Gt0_d_base[imy + N*L*jmy]);

					// Raman B1g

					const double delta_tau_j_ipx = (l == 0 && j == ipx ? 1 : 0);
					const double delta_tau_j_imx = (l == 0 && j == imx ? 1 : 0);
					const double delta_tau_j_ipy = (l == 0 && j == ipy ? 1 : 0);
					const double delta_tau_j_imy = (l == 0 && j == imy ? 1 : 0);

					meas_data->ram_b1g[i + Ncell*j + offset] += sign * 0.0625 * (
						+ ((Gtt_u_base[i + N*ipx] + Gtt_d_base[i + N*ipx])*(G00_u_base[j + N*jpx] + G00_d_base[j + N*jpx]) - Gt0_u_base[i + N*L*jpx]*(G0t_u_base[j + N*ipx] - delta_tau_j_ipx) - Gt0_d_base[i + N*L*jpx]*(G0t_d_base[j + N*ipx] - delta_tau_j_ipx))
						+ ((Gtt_u_base[i + N*ipx] + Gtt_d_base[i + N*ipx])*(G00_u_base[j + N*jmx] + G00_d_base[j + N*jmx]) - Gt0_u_base[i + N*L*jmx]*(G0t_u_base[j + N*ipx] - delta_tau_j_ipx) - Gt0_d_base[i + N*L*jmx]*(G0t_d_base[j + N*ipx] - delta_tau_j_ipx))
						- ((Gtt_u_base[i + N*ipx] + Gtt_d_base[i + N*ipx])*(G00_u_base[j + N*jpy] + G00_d_base[j + N*jpy]) - Gt0_u_base[i + N*L*jpy]*(G0t_u_base[j + N*ipx] - delta_tau_j_ipx) - Gt0_d_base[i + N*L*jpy]*(G0t_d_base[j + N*ipx] - delta_tau_j_ipx))
						- ((Gtt_u_base[i + N*ipx] + Gtt_d_base[i + N*ipx])*(G00_u_base[j + N*jmy] + G00_d_base[j + N*jmy]) - Gt0_u_base[i + N*L*jmy]*(G0t_u_base[j + N*ipx] - delta_tau_j_ipx) - Gt0_d_base[i + N*L*jmy]*(G0t_d_base[j + N*ipx] - delta_tau_j_ipx))

						+ ((Gtt_u_base[i + N*imx] + Gtt_d_base[i + N*imx])*(G00_u_base[j + N*jpx] + G00_d_base[j + N*jpx]) - Gt0_u_base[i + N*L*jpx]*(G0t_u_base[j + N*imx] - delta_tau_j_imx) - Gt0_d_base[i + N*L*jpx]*(G0t_d_base[j + N*imx] - delta_tau_j_imx))
						+ ((Gtt_u_base[i + N*imx] + Gtt_d_base[i + N*imx])*(G00_u_base[j + N*jmx] + G00_d_base[j + N*jmx]) - Gt0_u_base[i + N*L*jmx]*(G0t_u_base[j + N*imx] - delta_tau_j_imx) - Gt0_d_base[i + N*L*jmx]*(G0t_d_base[j + N*imx] - delta_tau_j_imx))
						- ((Gtt_u_base[i + N*imx] + Gtt_d_base[i + N*imx])*(G00_u_base[j + N*jpy] + G00_d_base[j + N*jpy]) - Gt0_u_base[i + N*L*jpy]*(G0t_u_base[j + N*imx] - delta_tau_j_imx) - Gt0_d_base[i + N*L*jpy]*(G0t_d_base[j + N*imx] - delta_tau_j_imx))
						- ((Gtt_u_base[i + N*imx] + Gtt_d_base[i + N*imx])*(G00_u_base[j + N*jmy] + G00_d_base[j + N*jmy]) - Gt0_u_base[i + N*L*jmy]*(G0t_u_base[j + N*imx] - delta_tau_j_imx) - Gt0_d_base[i + N*L*jmy]*(G0t_d_base[j + N*imx] - delta_tau_j_imx))

						- ((Gtt_u_base[i + N*ipy] + Gtt_d_base[i + N*ipy])*(G00_u_base[j + N*jpx] + G00_d_base[j + N*jpx]) - Gt0_u_base[i + N*L*jpx]*(G0t_u_base[j + N*ipy] - delta_tau_j_ipy) - Gt0_d_base[i + N*L*jpx]*(G0t_d_base[j + N*ipy] - delta_tau_j_ipy))
						- ((Gtt_u_base[i + N*ipy] + Gtt_d_base[i + N*ipy])*(G00_u_base[j + N*jmx] + G00_d_base[j + N*jmx]) - Gt0_u_base[i + N*L*jmx]*(G0t_u_base[j + N*ipy] - delta_tau_j_ipy) - Gt0_d_base[i + N*L*jmx]*(G0t_d_base[j + N*ipy] - delta_tau_j_ipy))
						+ ((Gtt_u_base[i + N*ipy] + Gtt_d_base[i + N*ipy])*(G00_u_base[j + N*jpy] + G00_d_base[j + N*jpy]) - Gt0_u_base[i + N*L*jpy]*(G0t_u_base[j + N*ipy] - delta_tau_j_ipy) - Gt0_d_base[i + N*L*jpy]*(G0t_d_base[j + N*ipy] - delta_tau_j_ipy))
						+ ((Gtt_u_base[i + N*ipy] + Gtt_d_base[i + N*ipy])*(G00_u_base[j + N*jmy] + G00_d_base[j + N*jmy]) - Gt0_u_base[i + N*L*jmy]*(G0t_u_base[j + N*ipy] - delta_tau_j_ipy) - Gt0_d_base[i + N*L*jmy]*(G0t_d_base[j + N*ipy] - delta_tau_j_ipy))

						- ((Gtt_u_base[i + N*imy] + Gtt_d_base[i + N*imy])*(G00_u_base[j + N*jpx] + G00_d_base[j + N*jpx]) - Gt0_u_base[i + N*L*jpx]*(G0t_u_base[j + N*imy] - delta_tau_j_imy) - Gt0_d_base[i + N*L*jpx]*(G0t_d_base[j + N*imy] - delta_tau_j_imy))
						- ((Gtt_u_base[i + N*imy] + Gtt_d_base[i + N*imy])*(G00_u_base[j + N*jmx] + G00_d_base[j + N*jmx]) - Gt0_u_base[i + N*L*jmx]*(G0t_u_base[j + N*imy] - delta_tau_j_imy) - Gt0_d_base[i + N*L*jmx]*(G0t_d_base[j + N*imy] - delta_tau_j_imy))
						+ ((Gtt_u_base[i + N*imy] + Gtt_d_base[i + N*imy])*(G00_u_base[j + N*jpy] + G00_d_base[j + N*jpy]) - Gt0_u_base[i + N*L*jpy]*(G0t_u_base[j + N*imy] - delta_tau_j_imy) - Gt0_d_base[i + N*L*jpy]*(G0t_d_base[j + N*imy] - delta_tau_j_imy))
						+ ((Gtt_u_base[i + N*imy] + Gtt_d_base[i + N*imy])*(G00_u_base[j + N*jmy] + G00_d_base[j + N*jmy]) - Gt0_u_base[i + N*L*jmy]*(G0t_u_base[j + N*imy] - delta_tau_j_imy) - Gt0_d_base[i + N*L*jmy]*(G0t_d_base[j + N*imy] - delta_tau_j_imy))
					);

					// Raman B2g

					const double delta_tau_j_ipxpy = (l == 0 && j == ipxpy ? 1 : 0);
					const double delta_tau_j_imxpy = (l == 0 && j == imxpy ? 1 : 0);
					const double delta_tau_j_ipxmy = (l == 0 && j == ipxmy ? 1 : 0);
					const double delta_tau_j_imxmy = (l == 0 && j == imxmy ? 1 : 0);

					meas_data->ram_b2g[i + Ncell*j + offset] += sign * 0.0625 * (
						+ ((Gtt_u_base[i + N*ipxpy] + Gtt_d_base[i + N*ipxpy])*(G00_u_base[j + N*jpxpy] + G00_d_base[j + N*jpxpy]) - Gt0_u_base[i + N*L*jpxpy]*(G0t_u_base[j + N*ipxpy] - delta_tau_j_ipxpy) - Gt0_d_base[i + N*L*jpxpy]*(G0t_d_base[j + N*ipxpy] - delta_tau_j_ipxpy))
						- ((Gtt_u_base[i + N*ipxpy] + Gtt_d_base[i + N*ipxpy])*(G00_u_base[j + N*jmxpy] + G00_d_base[j + N*jmxpy]) - Gt0_u_base[i + N*L*jmxpy]*(G0t_u_base[j + N*ipxpy] - delta_tau_j_ipxpy) - Gt0_d_base[i + N*L*jmxpy]*(G0t_d_base[j + N*ipxpy] - delta_tau_j_ipxpy))
						- ((Gtt_u_base[i + N*ipxpy] + Gtt_d_base[i + N*ipxpy])*(G00_u_base[j + N*jpxmy] + G00_d_base[j + N*jpxmy]) - Gt0_u_base[i + N*L*jpxmy]*(G0t_u_base[j + N*ipxpy] - delta_tau_j_ipxpy) - Gt0_d_base[i + N*L*jpxmy]*(G0t_d_base[j + N*ipxpy] - delta_tau_j_ipxpy))
						+ ((Gtt_u_base[i + N*ipxpy] + Gtt_d_base[i + N*ipxpy])*(G00_u_base[j + N*jmxmy] + G00_d_base[j + N*jmxmy]) - Gt0_u_base[i + N*L*jmxmy]*(G0t_u_base[j + N*ipxpy] - delta_tau_j_ipxpy) - Gt0_d_base[i + N*L*jmxmy]*(G0t_d_base[j + N*ipxpy] - delta_tau_j_ipxpy))

						- ((Gtt_u_base[i + N*imxpy] + Gtt_d_base[i + N*imxpy])*(G00_u_base[j + N*jpxpy] + G00_d_base[j + N*jpxpy]) - Gt0_u_base[i + N*L*jpxpy]*(G0t_u_base[j + N*imxpy] - delta_tau_j_imxpy) - Gt0_d_base[i + N*L*jpxpy]*(G0t_d_base[j + N*imxpy] - delta_tau_j_imxpy))
						+ ((Gtt_u_base[i + N*imxpy] + Gtt_d_base[i + N*imxpy])*(G00_u_base[j + N*jmxpy] + G00_d_base[j + N*jmxpy]) - Gt0_u_base[i + N*L*jmxpy]*(G0t_u_base[j + N*imxpy] - delta_tau_j_imxpy) - Gt0_d_base[i + N*L*jmxpy]*(G0t_d_base[j + N*imxpy] - delta_tau_j_imxpy))
						+ ((Gtt_u_base[i + N*imxpy] + Gtt_d_base[i + N*imxpy])*(G00_u_base[j + N*jpxmy] + G00_d_base[j + N*jpxmy]) - Gt0_u_base[i + N*L*jpxmy]*(G0t_u_base[j + N*imxpy] - delta_tau_j_imxpy) - Gt0_d_base[i + N*L*jpxmy]*(G0t_d_base[j + N*imxpy] - delta_tau_j_imxpy))
						- ((Gtt_u_base[i + N*imxpy] + Gtt_d_base[i + N*imxpy])*(G00_u_base[j + N*jmxmy] + G00_d_base[j + N*jmxmy]) - Gt0_u_base[i + N*L*jmxmy]*(G0t_u_base[j + N*imxpy] - delta_tau_j_imxpy) - Gt0_d_base[i + N*L*jmxmy]*(G0t_d_base[j + N*imxpy] - delta_tau_j_imxpy))

						- ((Gtt_u_base[i + N*ipxmy] + Gtt_d_base[i + N*ipxmy])*(G00_u_base[j + N*jpxpy] + G00_d_base[j + N*jpxpy]) - Gt0_u_base[i + N*L*jpxpy]*(G0t_u_base[j + N*ipxmy] - delta_tau_j_ipxmy) - Gt0_d_base[i + N*L*jpxpy]*(G0t_d_base[j + N*ipxmy] - delta_tau_j_ipxmy))
						+ ((Gtt_u_base[i + N*ipxmy] + Gtt_d_base[i + N*ipxmy])*(G00_u_base[j + N*jmxpy] + G00_d_base[j + N*jmxpy]) - Gt0_u_base[i + N*L*jmxpy]*(G0t_u_base[j + N*ipxmy] - delta_tau_j_ipxmy) - Gt0_d_base[i + N*L*jmxpy]*(G0t_d_base[j + N*ipxmy] - delta_tau_j_ipxmy))
						+ ((Gtt_u_base[i + N*ipxmy] + Gtt_d_base[i + N*ipxmy])*(G00_u_base[j + N*jpxmy] + G00_d_base[j + N*jpxmy]) - Gt0_u_base[i + N*L*jpxmy]*(G0t_u_base[j + N*ipxmy] - delta_tau_j_ipxmy) - Gt0_d_base[i + N*L*jpxmy]*(G0t_d_base[j + N*ipxmy] - delta_tau_j_ipxmy))
						- ((Gtt_u_base[i + N*ipxmy] + Gtt_d_base[i + N*ipxmy])*(G00_u_base[j + N*jmxmy] + G00_d_base[j + N*jmxmy]) - Gt0_u_base[i + N*L*jmxmy]*(G0t_u_base[j + N*ipxmy] - delta_tau_j_ipxmy) - Gt0_d_base[i + N*L*jmxmy]*(G0t_d_base[j + N*ipxmy] - delta_tau_j_ipxmy))

						+ ((Gtt_u_base[i + N*imxmy] + Gtt_d_base[i + N*imxmy])*(G00_u_base[j + N*jpxpy] + G00_d_base[j + N*jpxpy]) - Gt0_u_base[i + N*L*jpxpy]*(G0t_u_base[j + N*imxmy] - delta_tau_j_imxmy) - Gt0_d_base[i + N*L*jpxpy]*(G0t_d_base[j + N*imxmy] - delta_tau_j_imxmy))
						- ((Gtt_u_base[i + N*imxmy] + Gtt_d_base[i + N*imxmy])*(G00_u_base[j + N*jmxpy] + G00_d_base[j + N*jmxpy]) - Gt0_u_base[i + N*L*jmxpy]*(G0t_u_base[j + N*imxmy] - delta_tau_j_imxmy) - Gt0_d_base[i + N*L*jmxpy]*(G0t_d_base[j + N*imxmy] - delta_tau_j_imxmy))
						- ((Gtt_u_base[i + N*imxmy] + Gtt_d_base[i + N*imxmy])*(G00_u_base[j + N*jpxmy] + G00_d_base[j + N*jpxmy]) - Gt0_u_base[i + N*L*jpxmy]*(G0t_u_base[j + N*imxmy] - delta_tau_j_imxmy) - Gt0_d_base[i + N*L*jpxmy]*(G0t_d_base[j + N*imxmy] - delta_tau_j_imxmy))
						+ ((Gtt_u_base[i + N*imxmy] + Gtt_d_base[i + N*imxmy])*(G00_u_base[j + N*jmxmy] + G00_d_base[j + N*jmxmy]) - Gt0_u_base[i + N*L*jmxmy]*(G0t_u_base[j + N*imxmy] - delta_tau_j_imxmy) - Gt0_d_base[i + N*L*jmxmy]*(G0t_d_base[j + N*imxmy] - delta_tau_j_imxmy))
					);
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
/// \brief Normalize unequal time measurement data (divide by number of samples)
///
void NormalizeUnequalTimeMeasurementData(measurement_data_unequal_time_t *meas_data)
{
	// total number of orbitals and cells
	const int Norb  = meas_data->Norb;
	const int Ncell = meas_data->Ncell;

	const int N = Ncell * Norb;

	// total number of time steps
	const int L = meas_data->L;

	// normalization factor
	const double nfac = 1.0 / meas_data->nsampl;

	// divide by nsampl to get <s A>
	cblas_dscal(N*L*N, nfac, meas_data->Gtau0_u, 1);
	cblas_dscal(N*N*L, nfac, meas_data->G0tau_u, 1);
	cblas_dscal(N*N*L, nfac, meas_data->Geqlt_u, 1);
	cblas_dscal(N*L*N, nfac, meas_data->Gtau0_d, 1);
	cblas_dscal(N*N*L, nfac, meas_data->G0tau_d, 1);
	cblas_dscal(N*N*L, nfac, meas_data->Geqlt_d, 1);
	cblas_dscal(N*N*L, nfac, meas_data->nn_corr, 1);
	cblas_dscal(N*N*L, nfac, meas_data->zz_corr, 1);
	cblas_dscal(N*N*L, nfac, meas_data->xx_corr, 1);
	cblas_dscal(Ncell*Ncell*Norb*L, nfac, meas_data->sc_c_sw, 1);
	cblas_dscal(Ncell*Ncell*Norb*L, nfac, meas_data->sc_c_dw, 1);
	cblas_dscal(Ncell*Ncell*Norb*L, nfac, meas_data->sc_c_sx, 1);
	cblas_dscal(Ncell*Ncell*Norb*L, nfac, meas_data->ram_b1g, 1);
	cblas_dscal(Ncell*Ncell*Norb*L, nfac, meas_data->ram_b2g, 1);

	// calculate average sign
	meas_data->sign *= nfac;
}


//________________________________________________________________________________________________________________________
///
/// \brief Load unequal time measurement data structure from disk
///
void LoadUnequalTimeMeasurementData(const char *fnbase, const measurement_data_unequal_time_t *meas_data)
{
	const int Norb  = meas_data->Norb;
	const int Ncell = meas_data->Ncell;
	const int N = Ncell * Norb;
	const int L = meas_data->L;

	char path[1024];
	sprintf(path, "%s_uneqlt_sign.dat",    fnbase); ReadData(path, (void *)&meas_data->sign,   sizeof(double), 1);
	sprintf(path, "%s_uneqlt_nsampl.dat",  fnbase); ReadData(path, (void *)&meas_data->nsampl, sizeof(int), 1);

	sprintf(path, "%s_uneqlt_Gtau0_u.dat", fnbase); ReadData(path, meas_data->Gtau0_u, sizeof(double), N*L*N);
	sprintf(path, "%s_uneqlt_G0tau_u.dat", fnbase); ReadData(path, meas_data->G0tau_u, sizeof(double), N*N*L);
	sprintf(path, "%s_uneqlt_Geqlt_u.dat", fnbase); ReadData(path, meas_data->Geqlt_u, sizeof(double), N*N*L);
	sprintf(path, "%s_uneqlt_Gtau0_d.dat", fnbase); ReadData(path, meas_data->Gtau0_d, sizeof(double), N*L*N);
	sprintf(path, "%s_uneqlt_G0tau_d.dat", fnbase); ReadData(path, meas_data->G0tau_d, sizeof(double), N*N*L);
	sprintf(path, "%s_uneqlt_Geqlt_d.dat", fnbase); ReadData(path, meas_data->Geqlt_d, sizeof(double), N*N*L);

	sprintf(path, "%s_uneqlt_nn_corr.dat", fnbase); ReadData(path, meas_data->nn_corr, sizeof(double), N*N*L);
	sprintf(path, "%s_uneqlt_zz_corr.dat", fnbase); ReadData(path, meas_data->zz_corr, sizeof(double), N*N*L);
	sprintf(path, "%s_uneqlt_xx_corr.dat", fnbase); ReadData(path, meas_data->xx_corr, sizeof(double), N*N*L);

	sprintf(path, "%s_uneqlt_sc_c_sw.dat", fnbase); ReadData(path, meas_data->sc_c_sw, sizeof(double), Ncell*Ncell*Norb*L);
	sprintf(path, "%s_uneqlt_sc_c_dw.dat", fnbase); ReadData(path, meas_data->sc_c_dw, sizeof(double), Ncell*Ncell*Norb*L);
	sprintf(path, "%s_uneqlt_sc_c_sx.dat", fnbase); ReadData(path, meas_data->sc_c_sx, sizeof(double), Ncell*Ncell*Norb*L);

	sprintf(path, "%s_uneqlt_ram_b1g.dat", fnbase); ReadData(path, meas_data->ram_b1g, sizeof(double), Ncell*Ncell*Norb*L);
	sprintf(path, "%s_uneqlt_ram_b2g.dat", fnbase); ReadData(path, meas_data->ram_b2g, sizeof(double), Ncell*Ncell*Norb*L);
}


//________________________________________________________________________________________________________________________
///
/// \brief Save unequal time measurement data structure to disk
///
void SaveUnequalTimeMeasurementData(const char *fnbase, const measurement_data_unequal_time_t *meas_data)
{
	const int Norb  = meas_data->Norb;
	const int Ncell = meas_data->Ncell;
	const int N = Ncell * Norb;
	const int L = meas_data->L;

	char path[1024];
	sprintf(path, "%s_uneqlt_sign.dat",    fnbase); WriteData(path, &meas_data->sign,   sizeof(double), 1, false);
	sprintf(path, "%s_uneqlt_nsampl.dat",  fnbase); WriteData(path, &meas_data->nsampl, sizeof(int),    1, false);

	sprintf(path, "%s_uneqlt_Gtau0_u.dat", fnbase); WriteData(path, meas_data->Gtau0_u, sizeof(double), N*L*N, false);
	sprintf(path, "%s_uneqlt_G0tau_u.dat", fnbase); WriteData(path, meas_data->G0tau_u, sizeof(double), N*N*L, false);
	sprintf(path, "%s_uneqlt_Geqlt_u.dat", fnbase); WriteData(path, meas_data->Geqlt_u, sizeof(double), N*N*L, false);
	sprintf(path, "%s_uneqlt_Gtau0_d.dat", fnbase); WriteData(path, meas_data->Gtau0_d, sizeof(double), N*L*N, false);
	sprintf(path, "%s_uneqlt_G0tau_d.dat", fnbase); WriteData(path, meas_data->G0tau_d, sizeof(double), N*N*L, false);
	sprintf(path, "%s_uneqlt_Geqlt_d.dat", fnbase); WriteData(path, meas_data->Geqlt_d, sizeof(double), N*N*L, false);

	sprintf(path, "%s_uneqlt_nn_corr.dat", fnbase); WriteData(path, meas_data->nn_corr, sizeof(double), N*N*L, false);
	sprintf(path, "%s_uneqlt_zz_corr.dat", fnbase); WriteData(path, meas_data->zz_corr, sizeof(double), N*N*L, false);
	sprintf(path, "%s_uneqlt_xx_corr.dat", fnbase); WriteData(path, meas_data->xx_corr, sizeof(double), N*N*L, false);

	sprintf(path, "%s_uneqlt_sc_c_sw.dat", fnbase); WriteData(path, meas_data->sc_c_sw, sizeof(double), Ncell*Ncell*Norb*L, false);
	sprintf(path, "%s_uneqlt_sc_c_dw.dat", fnbase); WriteData(path, meas_data->sc_c_dw, sizeof(double), Ncell*Ncell*Norb*L, false);
	sprintf(path, "%s_uneqlt_sc_c_sx.dat", fnbase); WriteData(path, meas_data->sc_c_sx, sizeof(double), Ncell*Ncell*Norb*L, false);

	sprintf(path, "%s_uneqlt_ram_b1g.dat", fnbase); WriteData(path, meas_data->ram_b1g, sizeof(double), Ncell*Ncell*Norb*L, false);
	sprintf(path, "%s_uneqlt_ram_b2g.dat", fnbase); WriteData(path, meas_data->ram_b2g, sizeof(double), Ncell*Ncell*Norb*L, false);
}


//________________________________________________________________________________________________________________________
///
/// \brief Allocate and initialize phonon measurement data structure
///
void AllocatePhononData(const int Norb, const int Nx, const int Ny, const int pbc_shift, const int L, const int max_nsampl, measurement_data_phonon_t *restrict meas_data)
{
	// avoid "not referenced" compiler warning
	(void)pbc_shift;

	meas_data->Norb  = Norb;
	meas_data->Ncell = Nx * Ny;
	meas_data->L     = L;

	// no samples collected so far
	meas_data->nsampl = 0;
	meas_data->max_nsampl = max_nsampl;

	meas_data->sign = 0;

	meas_data->X_iteration = (double *)MKL_calloc(Norb * max_nsampl, sizeof(double), MEM_DATA_ALIGN);
	meas_data->X_avg       = (double *)MKL_calloc(Norb, sizeof(double), MEM_DATA_ALIGN);
	meas_data->X_sqr       = (double *)MKL_calloc(Norb, sizeof(double), MEM_DATA_ALIGN);
}


//________________________________________________________________________________________________________________________
///
/// \brief Delete phonon measurement data structure (free memory)
///
void DeletePhononData(measurement_data_phonon_t *restrict meas_data)
{
	MKL_free(meas_data->X_sqr);
	MKL_free(meas_data->X_avg);
	MKL_free(meas_data->X_iteration);
}


//________________________________________________________________________________________________________________________
///
/// \brief Accumulate phonon "measurement" data
///
void AccumulatePhononData(const greens_func_t *restrict Gu, const greens_func_t *restrict Gd, const double *restrict X, measurement_data_phonon_t *restrict meas_data)
{
	// lattice dimensions
	const int Norb  = meas_data->Norb;
	const int Ncell = meas_data->Ncell;
	const int L = meas_data->L;
	const int N = Ncell * Norb;

	// product of the determinant signs of the Green's function matrices
	const double sign = (double)(Gu->sgndet * Gd->sgndet);
	assert(sign != 0);

	// sign and normalization factor
	const double signfac = sign / Ncell / L;

	int l;
	for (l = 0; l < L; l++)		// for all discrete time differences...
	{
		int o;
		for (o = 0; o < Norb; o++)
		{
			int i;
			for (i = 0; i < Ncell; i++)
			{
				meas_data->X_avg[o] += signfac * X[i + o*Ncell + l*N];
				meas_data->X_sqr[o] += signfac * square(X[i + o*Ncell + l*N]);
				meas_data->X_iteration[o + meas_data->nsampl * Norb] += X[i + o*Ncell + l*N] / (Ncell * L);
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
/// \brief Normalize phonon measurement data (divide by number of samples)
///
void NormalizePhononData(measurement_data_phonon_t *meas_data)
{
	// total number of orbitals and cells
	const int Norb  = meas_data->Norb;

	// normalization factor
	const double nfac = 1.0 / meas_data->nsampl;

	// divide by nsampl to get <s A>
	cblas_dscal(Norb, nfac, meas_data->X_avg, 1);
	cblas_dscal(Norb, nfac, meas_data->X_sqr, 1);

	// calculate average sign
	meas_data->sign *= nfac;
}


//________________________________________________________________________________________________________________________
///
/// \brief Print a basic phonon measurement summary
///
void PrintPhononData(const measurement_data_phonon_t *meas_data)
{
	duprintf("_______________________________________________________________________________\n");
	duprintf("Summary of phonon data\n\n");
	duprintf("                    average sign: %g\n", meas_data->sign);
	int o;
	for (o = 0; o < meas_data->Norb; o++)
	{
		duprintf("\nResults for orbital %d\n", o);
		duprintf("                       average X: %g\n", meas_data->X_avg[o]);
		duprintf("                     average X^2: %g\n", meas_data->X_sqr[o]);
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief Save phonon measurement data structure to disk
///
void SavePhononData(const char *fnbase, const measurement_data_phonon_t *meas_data)
{
	const int Norb       = meas_data->Norb;
	const int max_nsampl = meas_data->max_nsampl;

	char path[1024];
	sprintf(path, "%s_phonon_sign.dat",        fnbase); WriteData(path, &meas_data->sign,       sizeof(double), 1, false);
	sprintf(path, "%s_phonon_nsampl.dat",      fnbase); WriteData(path, &meas_data->nsampl,     sizeof(int),    1, false);
	sprintf(path, "%s_phonon_X_iteration.dat", fnbase); WriteData(path, meas_data->X_iteration, sizeof(double), Norb * max_nsampl, false);
	sprintf(path, "%s_phonon_X_avg.dat",       fnbase); WriteData(path, meas_data->X_avg,       sizeof(double), Norb, false);
	sprintf(path, "%s_phonon_X_sq.dat",        fnbase); WriteData(path, meas_data->X_sqr,       sizeof(double), Norb, false);
}
