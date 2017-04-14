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
static void ConstructLatticeCoordinateSumMap(const int Nx, const int Ny, const int pbc_shift, int *restrict latt_sum_map)
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
					const int u = (i + k + (j + l < Ny ? 0 : pbc_shift)) % Nx;
					const int v = (j + l)                                % Ny;

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

			latt_xp1_map[ij] = i_next + Nx*j;           // index of (i+1,j) lattice site
			latt_xm1_map[ij] = i_prev + Nx*j;           // index of (i-1,j) lattice site
			latt_yp1_map[ij] = i      + Nx*j_next;      // index of (i,j+1) lattice site
			latt_ym1_map[ij] = i      + Nx*j_prev;      // index of (i,j-1) lattice site
		}
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief Allocate and initialize measurement data structure
///
void AllocateMeasurementData(const int Norb, const int Nx, const int Ny, const int pbc_shift, measurement_data_t *restrict meas_data)
{
	const int Ncell = Nx * Ny;
	const int N     = Ncell * Norb;

	meas_data->Norb  = Norb;
	meas_data->Ncell = Ncell;

	// no samples collected so far
	meas_data->nsampl = 0;

	meas_data->sign = 0;

	meas_data->density_u = (double *)MKL_calloc(Norb, sizeof(double), MEM_DATA_ALIGN);
	meas_data->density_d = (double *)MKL_calloc(Norb, sizeof(double), MEM_DATA_ALIGN);
	meas_data->doubleocc = (double *)MKL_calloc(Norb, sizeof(double), MEM_DATA_ALIGN);

	// green's functions
	meas_data->grfun_u = (double *)MKL_calloc(N*N, sizeof(double), MEM_DATA_ALIGN);
	meas_data->grfun_d = (double *)MKL_calloc(N*N, sizeof(double), MEM_DATA_ALIGN);

	// density correlation data
	// first index for spatial offset, second and third indices for orbital index
	meas_data->uu_corr = (double *)MKL_calloc(Ncell*Norb*Norb, sizeof(double), MEM_DATA_ALIGN);
	meas_data->dd_corr = (double *)MKL_calloc(Ncell*Norb*Norb, sizeof(double), MEM_DATA_ALIGN);
	meas_data->ud_corr = (double *)MKL_calloc(Ncell*Norb*Norb, sizeof(double), MEM_DATA_ALIGN);
	meas_data->ff_corr = (double *)MKL_calloc(Ncell*Norb*Norb, sizeof(double), MEM_DATA_ALIGN);

	// spin correlation data
	meas_data->zz_corr = (double *)MKL_calloc(Ncell*Norb*Norb, sizeof(double), MEM_DATA_ALIGN);
	meas_data->xx_corr = (double *)MKL_calloc(Ncell*Norb*Norb, sizeof(double), MEM_DATA_ALIGN);

	// superconducting susceptibilities
	meas_data->sc_c_sw = (double *)MKL_calloc(N, sizeof(double), MEM_DATA_ALIGN);
	meas_data->sc_c_dw = (double *)MKL_calloc(N, sizeof(double), MEM_DATA_ALIGN);
	meas_data->sc_c_sx = (double *)MKL_calloc(N, sizeof(double), MEM_DATA_ALIGN);

	// construct lattice coordinate sum map
	meas_data->latt_sum_map = (int *)MKL_malloc(Ncell*Ncell * sizeof(int), MEM_DATA_ALIGN);
	ConstructLatticeCoordinateSumMap(Nx, Ny, pbc_shift, meas_data->latt_sum_map);

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

	MKL_free(meas_data->latt_sum_map);

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
	int i, k;   // spatial indices
	int o, p;   // orbital indices

	// lattice dimensions
	const int Norb  = meas_data->Norb;
	const int Ncell = meas_data->Ncell;
	const int N = Ncell * Norb;

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
			const int io = i + Ncell * o;
			nu += (1 - Gu->mat[io + N*io]);
			nd += (1 - Gd->mat[io + N*io]);
			oc += (1 - Gu->mat[io + N*io])*(1 - Gd->mat[io + N*io]);
		}
		// mean value of density and double occupancy for current orbital
		meas_data->density_u[o] += signfac*nu;
		meas_data->density_d[o] += signfac*nd;
		meas_data->doubleocc[o] += signfac*oc;
	}

	// Green's functions
	cblas_daxpy(N*N, sign, Gu->mat, 1, meas_data->grfun_u, 1);
	cblas_daxpy(N*N, sign, Gd->mat, 1, meas_data->grfun_d, 1);

	// density and spin correlations

	for (o = 0; o < Norb; o++)
	{
		for (p = 0; p < Norb; p++)
		{
			// starting index for correlation accumulators
			const int offset = Ncell * (o + Norb*p);

			for (i = 0; i < Ncell; i++)
			{
				const int io = i + Ncell * o;
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

					if (k == 0 && o == p)   // special case: same lattice site and orbital
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

	// superconducting susceptibilities (separately for each orbital type)

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

			for (k = 0; k < Ncell; k++)
			{
				const int j = meas_data->latt_sum_map[k + Ncell*i];

				// nearest neighbors of lattice site j
				const int jpx = meas_data->latt_xp1_map[j];
				const int jmx = meas_data->latt_xm1_map[j];
				const int jpy = meas_data->latt_yp1_map[j];
				const int jmy = meas_data->latt_ym1_map[j];

				const double Gu_ij = Gu_orb[i + N*j];
				const double Gd_ij = Gd_orb[i + N*j];

				// s-wave
				meas_data->sc_c_sw[k + Ncell*o] += signfac*(Gu_ij * Gd_ij);

				// d-wave: all 16 combinations of nearest neighbors of i and j
				meas_data->sc_c_dw[k + Ncell*o] += signfac * 0.25 * Gu_ij*(
					+ Gd_orb[ipx + N*jpx] + Gd_orb[ipx + N*jmx] - Gd_orb[ipx + N*jpy] - Gd_orb[ipx + N*jmy]
					+ Gd_orb[imx + N*jpx] + Gd_orb[imx + N*jmx] - Gd_orb[imx + N*jpy] - Gd_orb[imx + N*jmy]
					- Gd_orb[ipy + N*jpx] - Gd_orb[ipy + N*jmx] + Gd_orb[ipy + N*jpy] + Gd_orb[ipy + N*jmy]
					- Gd_orb[imy + N*jpx] - Gd_orb[imy + N*jmx] + Gd_orb[imy + N*jpy] + Gd_orb[imy + N*jmy]);

				// extended s-wave: similar to d-wave, but without sign flip
				meas_data->sc_c_sx[k + Ncell*o] += signfac * 0.25 * Gu_ij*(
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
	cblas_dscal(Norb, nfac, meas_data->density_u, 1);
	cblas_dscal(Norb, nfac, meas_data->density_d, 1);
	cblas_dscal(Norb, nfac, meas_data->doubleocc, 1);
	cblas_dscal(N*N, nfac, meas_data->grfun_u, 1);
	cblas_dscal(N*N, nfac, meas_data->grfun_d, 1);
	const int m = Ncell*Norb*Norb;
	cblas_dscal(m, nfac, meas_data->uu_corr, 1);
	cblas_dscal(m, nfac, meas_data->dd_corr, 1);
	cblas_dscal(m, nfac, meas_data->ud_corr, 1);
	cblas_dscal(m, nfac, meas_data->ff_corr, 1);
	cblas_dscal(m, nfac, meas_data->zz_corr, 1);
	cblas_dscal(m, nfac, meas_data->xx_corr, 1);
	cblas_dscal(N, nfac, meas_data->sc_c_sw, 1);
	cblas_dscal(N, nfac, meas_data->sc_c_dw, 1);
	cblas_dscal(N, nfac, meas_data->sc_c_sx, 1);

	// calculate average sign
	meas_data->sign *= nfac;
}


//________________________________________________________________________________________________________________________
///
/// \brief Print a summary of basic measurements
///
void PrintMeasurementDataSummary(const measurement_data_t *meas_data)
{
	const double sgn = meas_data->sign;
	duprintf("_______________________________________________________________________________\n");
	duprintf("Summary of simulation results\n\n");
	duprintf("                    average sign: %g\n", sgn);
	double total_density = 0;
	int o;
	for (o = 0; o < meas_data->Norb; o++)
	{
		total_density += meas_data->density_u[o] + meas_data->density_d[o];
	}
	duprintf("           average total density: %g\n", total_density/sgn);
	for (o = 0; o < meas_data->Norb; o++)
	{
		duprintf("\nResults for orbital %d\n", o);
		duprintf("           average total density: %g\n", (meas_data->density_u[o] + meas_data->density_d[o])/sgn);
		duprintf("         average spin-up density: %g\n", meas_data->density_u[o]/sgn);
		duprintf("       average spin-down density: %g\n", meas_data->density_d[o]/sgn);
		duprintf("        average double occupancy: %g\n", meas_data->doubleocc[o]/sgn);
		duprintf("            average local moment: %g\n", (meas_data->density_u[o] + meas_data->density_d[o] - 2.0*meas_data->doubleocc[o])/sgn);
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
	const int N     = Ncell * Norb;

	char path[1024];
	sprintf(path, "%s_eqlt_sign.dat",      fnbase); ReadData(path, (void *)&meas_data->sign,   sizeof(double), 1);
	sprintf(path, "%s_eqlt_nsampl.dat",    fnbase); ReadData(path, (void *)&meas_data->nsampl, sizeof(int),    1);
	sprintf(path, "%s_eqlt_density_u.dat", fnbase); ReadData(path, meas_data->density_u, sizeof(double), Norb);
	sprintf(path, "%s_eqlt_density_d.dat", fnbase); ReadData(path, meas_data->density_d, sizeof(double), Norb);
	sprintf(path, "%s_eqlt_doubleocc.dat", fnbase); ReadData(path, meas_data->doubleocc, sizeof(double), Norb);
	sprintf(path, "%s_eqlt_grfun_u.dat",   fnbase); ReadData(path, meas_data->grfun_u,   sizeof(double), N*N);
	sprintf(path, "%s_eqlt_grfun_d.dat",   fnbase); ReadData(path, meas_data->grfun_d,   sizeof(double), N*N);
	sprintf(path, "%s_eqlt_uu_corr.dat",   fnbase); ReadData(path, meas_data->uu_corr,   sizeof(double), Ncell*Norb*Norb);
	sprintf(path, "%s_eqlt_dd_corr.dat",   fnbase); ReadData(path, meas_data->dd_corr,   sizeof(double), Ncell*Norb*Norb);
	sprintf(path, "%s_eqlt_ud_corr.dat",   fnbase); ReadData(path, meas_data->ud_corr,   sizeof(double), Ncell*Norb*Norb);
	sprintf(path, "%s_eqlt_ff_corr.dat",   fnbase); ReadData(path, meas_data->ff_corr,   sizeof(double), Ncell*Norb*Norb);
	sprintf(path, "%s_eqlt_zz_corr.dat",   fnbase); ReadData(path, meas_data->zz_corr,   sizeof(double), Ncell*Norb*Norb);
	sprintf(path, "%s_eqlt_xx_corr.dat",   fnbase); ReadData(path, meas_data->xx_corr,   sizeof(double), Ncell*Norb*Norb);
	sprintf(path, "%s_eqlt_sc_c_sw.dat",   fnbase); ReadData(path, meas_data->sc_c_sw,   sizeof(double), N);
	sprintf(path, "%s_eqlt_sc_c_dw.dat",   fnbase); ReadData(path, meas_data->sc_c_dw,   sizeof(double), N);
	sprintf(path, "%s_eqlt_sc_c_sx.dat",   fnbase); ReadData(path, meas_data->sc_c_sx,   sizeof(double), N);
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
	sprintf(path, "%s_eqlt_density_u.dat", fnbase); WriteData(path, meas_data->density_u, sizeof(double), Norb, false);
	sprintf(path, "%s_eqlt_density_d.dat", fnbase); WriteData(path, meas_data->density_d, sizeof(double), Norb, false);
	sprintf(path, "%s_eqlt_doubleocc.dat", fnbase); WriteData(path, meas_data->doubleocc, sizeof(double), Norb, false);
	sprintf(path, "%s_eqlt_grfun_u.dat",   fnbase); WriteData(path, meas_data->grfun_u,   sizeof(double), N*N, false);
	sprintf(path, "%s_eqlt_grfun_d.dat",   fnbase); WriteData(path, meas_data->grfun_d,   sizeof(double), N*N, false);
	sprintf(path, "%s_eqlt_uu_corr.dat",   fnbase); WriteData(path, meas_data->uu_corr,   sizeof(double), Ncell*Norb*Norb, false);
	sprintf(path, "%s_eqlt_dd_corr.dat",   fnbase); WriteData(path, meas_data->dd_corr,   sizeof(double), Ncell*Norb*Norb, false);
	sprintf(path, "%s_eqlt_ud_corr.dat",   fnbase); WriteData(path, meas_data->ud_corr,   sizeof(double), Ncell*Norb*Norb, false);
	sprintf(path, "%s_eqlt_ff_corr.dat",   fnbase); WriteData(path, meas_data->ff_corr,   sizeof(double), Ncell*Norb*Norb, false);
	sprintf(path, "%s_eqlt_zz_corr.dat",   fnbase); WriteData(path, meas_data->zz_corr,   sizeof(double), Ncell*Norb*Norb, false);
	sprintf(path, "%s_eqlt_xx_corr.dat",   fnbase); WriteData(path, meas_data->xx_corr,   sizeof(double), Ncell*Norb*Norb, false);
	sprintf(path, "%s_eqlt_sc_c_sw.dat",   fnbase); WriteData(path, meas_data->sc_c_sw,   sizeof(double), N, false);
	sprintf(path, "%s_eqlt_sc_c_dw.dat",   fnbase); WriteData(path, meas_data->sc_c_dw,   sizeof(double), N, false);
	sprintf(path, "%s_eqlt_sc_c_sx.dat",   fnbase); WriteData(path, meas_data->sc_c_sx,   sizeof(double), N, false);
}


//________________________________________________________________________________________________________________________
///
/// \brief Allocate and initialize unequal time measurement data structure
///
int AllocateUnequalTimeMeasurementData(const int Norb, const int Nx, const int Ny, const int pbc_shift, const int L, const int numBprod, measurement_data_unequal_time_t *restrict meas_data)
{
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
	meas_data->Gtau0_u = (double *)MKL_calloc(Ncell*Norb*Norb*L, sizeof(double), MEM_DATA_ALIGN);
	meas_data->G0tau_u = (double *)MKL_calloc(Ncell*Norb*Norb*L, sizeof(double), MEM_DATA_ALIGN);
	meas_data->Geqlt_u = (double *)MKL_calloc(Ncell*Norb*Norb*L, sizeof(double), MEM_DATA_ALIGN);
	meas_data->Gtau0_d = (double *)MKL_calloc(Ncell*Norb*Norb*L, sizeof(double), MEM_DATA_ALIGN);
	meas_data->G0tau_d = (double *)MKL_calloc(Ncell*Norb*Norb*L, sizeof(double), MEM_DATA_ALIGN);
	meas_data->Geqlt_d = (double *)MKL_calloc(Ncell*Norb*Norb*L, sizeof(double), MEM_DATA_ALIGN);

	// density and spin correlation data for all time differences
	meas_data->nn_corr = (double *)MKL_calloc(Ncell*Norb*Norb*L, sizeof(double), MEM_DATA_ALIGN);
	meas_data->zz_corr = (double *)MKL_calloc(Ncell*Norb*Norb*L, sizeof(double), MEM_DATA_ALIGN);
	meas_data->xx_corr = (double *)MKL_calloc(Ncell*Norb*Norb*L, sizeof(double), MEM_DATA_ALIGN);

	// superconducting susceptibilities
	meas_data->sc_c_sw = (double *)MKL_calloc(N*L, sizeof(double), MEM_DATA_ALIGN);
	meas_data->sc_c_dw = (double *)MKL_calloc(N*L, sizeof(double), MEM_DATA_ALIGN);
	meas_data->sc_c_sx = (double *)MKL_calloc(N*L, sizeof(double), MEM_DATA_ALIGN);

	// Raman
	meas_data->ram_b1g = (double *)MKL_calloc(N*L, sizeof(double), MEM_DATA_ALIGN);
	meas_data->ram_b2g = (double *)MKL_calloc(N*L, sizeof(double), MEM_DATA_ALIGN);

	// construct lattice coordinate sum map
	meas_data->latt_sum_map = (int *)MKL_malloc(Ncell*Ncell * sizeof(int), MEM_DATA_ALIGN);
	ConstructLatticeCoordinateSumMap(Nx, Ny, pbc_shift, meas_data->latt_sum_map);

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
	int i, k;
	int o, p;

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

	// compute unequal time spin-up and spin-down Green's functions
	#pragma omp parallel sections
	{
		#pragma omp section
		ComputeUnequalTimeGreensFunction(N, L, tsm_u, meas_data->Hu, curGtau0_u, curG0tau_u, curGeqlt_u);
		#pragma omp section
		ComputeUnequalTimeGreensFunction(N, L, tsm_d, meas_data->Hd, curGtau0_d, curG0tau_d, curGeqlt_d);
	}

	// accumulate density and spin correlation data, as well as (unequal-time) Green's functions

	// sign and normalization factor
	const double signfac = sign / Ncell;

	for (l = 0; l < L; l++)     // for all discrete time differences...
	{
		for (o = 0; o < Norb; o++)
		{
			for (p = 0; p < Norb; p++)
			{
				const int offset = Ncell * (o + Norb*p + Norb*Norb*l);

				for (i = 0; i < Ncell; i++)
				{
					const int io = i + o * Ncell;
					const double Gtt_u_ii = curGeqlt_u[io + N*(io + N*l)];
					const double Gtt_d_ii = curGeqlt_d[io + N*(io + N*l)];

					for (k = 0; k < Ncell; k++)
					{
						const int jp = meas_data->latt_sum_map[k + Ncell*i] + Ncell*p;
						meas_data->Geqlt_u[k + offset] += signfac*curGeqlt_u[io + N*(jp + N*l)];
						meas_data->Geqlt_d[k + offset] += signfac*curGeqlt_d[io + N*(jp + N*l)];
						meas_data->Gtau0_u[k + offset] += signfac*curGtau0_u[io + N*(l + L*jp)];
						meas_data->Gtau0_d[k + offset] += signfac*curGtau0_d[io + N*(l + L*jp)];
						meas_data->G0tau_u[k + offset] += signfac*curG0tau_u[io + N*(jp + N*l)];
						meas_data->G0tau_d[k + offset] += signfac*curG0tau_d[io + N*(jp + N*l)];

						// special case: same lattice site, orbital and time slice
						if (k == 0 && o == p && l == 0)
						{
							meas_data->nn_corr[0 + offset] += signfac*(4 - 3*(Gtt_u_ii + Gtt_d_ii) + 2*Gtt_u_ii*Gtt_d_ii);
							meas_data->zz_corr[0 + offset] += signfac*(Gtt_u_ii - 2*Gtt_u_ii*Gtt_d_ii + Gtt_d_ii);
							meas_data->xx_corr[0 + offset] += signfac*(Gtt_u_ii - 2*Gtt_u_ii*Gtt_d_ii + Gtt_d_ii);
						}
						else
						{
							const int jp = meas_data->latt_sum_map[k + Ncell*i] + Ncell*p;

							const double G00_u_jj = curGeqlt_u[jp + N*jp];
							const double G00_d_jj = curGeqlt_d[jp + N*jp];

							const double Gt0_u_ij = curGtau0_u[io + N*(l + L*jp)];
							const double Gt0_d_ij = curGtau0_d[io + N*(l + L*jp)];
							const double G0t_u_ji = curG0tau_u[jp + N*(io + N*l)];
							const double G0t_d_ji = curG0tau_d[jp + N*(io + N*l)];

							meas_data->nn_corr[k + offset] += signfac*((2 - Gtt_u_ii - Gtt_d_ii)*(2 - G00_u_jj - G00_d_jj) - (Gt0_u_ij*G0t_u_ji + Gt0_d_ij*G0t_d_ji));
							meas_data->zz_corr[k + offset] += signfac*((    Gtt_u_ii - Gtt_d_ii)*(    G00_u_jj - G00_d_jj) - (Gt0_u_ij*G0t_u_ji + Gt0_d_ij*G0t_d_ji));
							meas_data->xx_corr[k + offset] += signfac*(                                                    - (Gt0_u_ij*G0t_d_ji + Gt0_d_ij*G0t_u_ji));
						}
					}
				}
			}


			// superconducting susceptibilities and Raman correlation functions (separately for each orbital type)

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

			const int offset = Ncell*o + N*l;

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

				for (k = 0; k < Ncell; k++)
				{
					const int j = meas_data->latt_sum_map[k + Ncell*i];

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
					meas_data->sc_c_sw[k + offset] += signfac*(Gt0_u_ij * Gt0_d_ij);

					// d-wave: all 16 combinations of nearest neighbors of i and j
					meas_data->sc_c_dw[k + offset] += signfac * 0.25 * Gt0_u_ij*(
						+ Gt0_d_base[ipx + N*L*jpx] + Gt0_d_base[ipx + N*L*jmx] - Gt0_d_base[ipx + N*L*jpy] - Gt0_d_base[ipx + N*L*jmy]
						+ Gt0_d_base[imx + N*L*jpx] + Gt0_d_base[imx + N*L*jmx] - Gt0_d_base[imx + N*L*jpy] - Gt0_d_base[imx + N*L*jmy]
						- Gt0_d_base[ipy + N*L*jpx] - Gt0_d_base[ipy + N*L*jmx] + Gt0_d_base[ipy + N*L*jpy] + Gt0_d_base[ipy + N*L*jmy]
						- Gt0_d_base[imy + N*L*jpx] - Gt0_d_base[imy + N*L*jmx] + Gt0_d_base[imy + N*L*jpy] + Gt0_d_base[imy + N*L*jmy]);

					// extended s-wave: similar to d-wave, but without sign flip
					meas_data->sc_c_sx[k + offset] += signfac * 0.25 * Gt0_u_ij*(
						+ Gt0_d_base[ipx + N*L*jpx] + Gt0_d_base[ipx + N*L*jmx] + Gt0_d_base[ipx + N*L*jpy] + Gt0_d_base[ipx + N*L*jmy]
						+ Gt0_d_base[imx + N*L*jpx] + Gt0_d_base[imx + N*L*jmx] + Gt0_d_base[imx + N*L*jpy] + Gt0_d_base[imx + N*L*jmy]
						+ Gt0_d_base[ipy + N*L*jpx] + Gt0_d_base[ipy + N*L*jmx] + Gt0_d_base[ipy + N*L*jpy] + Gt0_d_base[ipy + N*L*jmy]
						+ Gt0_d_base[imy + N*L*jpx] + Gt0_d_base[imy + N*L*jmx] + Gt0_d_base[imy + N*L*jpy] + Gt0_d_base[imy + N*L*jmy]);

					// Raman B1g

					const double delta_tau_j_ipx = (l == 0 && j == ipx ? 1 : 0);
					const double delta_tau_j_imx = (l == 0 && j == imx ? 1 : 0);
					const double delta_tau_j_ipy = (l == 0 && j == ipy ? 1 : 0);
					const double delta_tau_j_imy = (l == 0 && j == imy ? 1 : 0);

					meas_data->ram_b1g[k + offset] += signfac * 0.0625 * (
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

					meas_data->ram_b2g[k + offset] += signfac * 0.0625 * (
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
	cblas_dscal(Ncell*Norb*Norb*L, nfac, meas_data->Gtau0_u, 1);
	cblas_dscal(Ncell*Norb*Norb*L, nfac, meas_data->G0tau_u, 1);
	cblas_dscal(Ncell*Norb*Norb*L, nfac, meas_data->Geqlt_u, 1);
	cblas_dscal(Ncell*Norb*Norb*L, nfac, meas_data->Gtau0_d, 1);
	cblas_dscal(Ncell*Norb*Norb*L, nfac, meas_data->G0tau_d, 1);
	cblas_dscal(Ncell*Norb*Norb*L, nfac, meas_data->Geqlt_d, 1);
	cblas_dscal(Ncell*Norb*Norb*L, nfac, meas_data->nn_corr, 1);
	cblas_dscal(Ncell*Norb*Norb*L, nfac, meas_data->zz_corr, 1);
	cblas_dscal(Ncell*Norb*Norb*L, nfac, meas_data->xx_corr, 1);
	cblas_dscal(N*L, nfac, meas_data->sc_c_sw, 1);
	cblas_dscal(N*L, nfac, meas_data->sc_c_dw, 1);
	cblas_dscal(N*L, nfac, meas_data->sc_c_sx, 1);
	cblas_dscal(N*L, nfac, meas_data->ram_b1g, 1);
	cblas_dscal(N*L, nfac, meas_data->ram_b2g, 1);

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

	sprintf(path, "%s_uneqlt_Gtau0_u.dat", fnbase); ReadData(path, meas_data->Gtau0_u, sizeof(double), Ncell*Norb*Norb*L);
	sprintf(path, "%s_uneqlt_G0tau_u.dat", fnbase); ReadData(path, meas_data->G0tau_u, sizeof(double), Ncell*Norb*Norb*L);
	sprintf(path, "%s_uneqlt_Geqlt_u.dat", fnbase); ReadData(path, meas_data->Geqlt_u, sizeof(double), Ncell*Norb*Norb*L);
	sprintf(path, "%s_uneqlt_Gtau0_d.dat", fnbase); ReadData(path, meas_data->Gtau0_d, sizeof(double), Ncell*Norb*Norb*L);
	sprintf(path, "%s_uneqlt_G0tau_d.dat", fnbase); ReadData(path, meas_data->G0tau_d, sizeof(double), Ncell*Norb*Norb*L);
	sprintf(path, "%s_uneqlt_Geqlt_d.dat", fnbase); ReadData(path, meas_data->Geqlt_d, sizeof(double), Ncell*Norb*Norb*L);

	sprintf(path, "%s_uneqlt_nn_corr.dat", fnbase); ReadData(path, meas_data->nn_corr, sizeof(double), Ncell*Norb*Norb*L);
	sprintf(path, "%s_uneqlt_zz_corr.dat", fnbase); ReadData(path, meas_data->zz_corr, sizeof(double), Ncell*Norb*Norb*L);
	sprintf(path, "%s_uneqlt_xx_corr.dat", fnbase); ReadData(path, meas_data->xx_corr, sizeof(double), Ncell*Norb*Norb*L);

	sprintf(path, "%s_uneqlt_sc_c_sw.dat", fnbase); ReadData(path, meas_data->sc_c_sw, sizeof(double), N*L);
	sprintf(path, "%s_uneqlt_sc_c_dw.dat", fnbase); ReadData(path, meas_data->sc_c_dw, sizeof(double), N*L);
	sprintf(path, "%s_uneqlt_sc_c_sx.dat", fnbase); ReadData(path, meas_data->sc_c_sx, sizeof(double), N*L);

	sprintf(path, "%s_uneqlt_ram_b1g.dat", fnbase); ReadData(path, meas_data->ram_b1g, sizeof(double), N*L);
	sprintf(path, "%s_uneqlt_ram_b2g.dat", fnbase); ReadData(path, meas_data->ram_b2g, sizeof(double), N*L);
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

	sprintf(path, "%s_uneqlt_Gtau0_u.dat", fnbase); WriteData(path, meas_data->Gtau0_u, sizeof(double), Ncell*Norb*Norb*L, false);
	sprintf(path, "%s_uneqlt_G0tau_u.dat", fnbase); WriteData(path, meas_data->G0tau_u, sizeof(double), Ncell*Norb*Norb*L, false);
	sprintf(path, "%s_uneqlt_Geqlt_u.dat", fnbase); WriteData(path, meas_data->Geqlt_u, sizeof(double), Ncell*Norb*Norb*L, false);
	sprintf(path, "%s_uneqlt_Gtau0_d.dat", fnbase); WriteData(path, meas_data->Gtau0_d, sizeof(double), Ncell*Norb*Norb*L, false);
	sprintf(path, "%s_uneqlt_G0tau_d.dat", fnbase); WriteData(path, meas_data->G0tau_d, sizeof(double), Ncell*Norb*Norb*L, false);
	sprintf(path, "%s_uneqlt_Geqlt_d.dat", fnbase); WriteData(path, meas_data->Geqlt_d, sizeof(double), Ncell*Norb*Norb*L, false);

	sprintf(path, "%s_uneqlt_nn_corr.dat", fnbase); WriteData(path, meas_data->nn_corr, sizeof(double), Ncell*Norb*Norb*L, false);
	sprintf(path, "%s_uneqlt_zz_corr.dat", fnbase); WriteData(path, meas_data->zz_corr, sizeof(double), Ncell*Norb*Norb*L, false);
	sprintf(path, "%s_uneqlt_xx_corr.dat", fnbase); WriteData(path, meas_data->xx_corr, sizeof(double), Ncell*Norb*Norb*L, false);

	sprintf(path, "%s_uneqlt_sc_c_sw.dat", fnbase); WriteData(path, meas_data->sc_c_sw, sizeof(double), N*L, false);
	sprintf(path, "%s_uneqlt_sc_c_dw.dat", fnbase); WriteData(path, meas_data->sc_c_dw, sizeof(double), N*L, false);
	sprintf(path, "%s_uneqlt_sc_c_sx.dat", fnbase); WriteData(path, meas_data->sc_c_sx, sizeof(double), N*L, false);

	sprintf(path, "%s_uneqlt_ram_b1g.dat", fnbase); WriteData(path, meas_data->ram_b1g, sizeof(double), N*L, false);
	sprintf(path, "%s_uneqlt_ram_b2g.dat", fnbase); WriteData(path, meas_data->ram_b2g, sizeof(double), N*L, false);
}


//________________________________________________________________________________________________________________________
///
/// \brief Allocate phonon measurement data
///
void AllocatePhononData(const int Norb, const int Nx, const int Ny, const int pbc_shift, const int L, const int max_nsampl, measurement_data_phonon_t *restrict meas_data)
{
	const int Ncell = Nx * Ny;

	meas_data->Norb  = Norb;
	meas_data->Ncell = Ncell;
	meas_data->L     = L;

	// no samples collected so far
	meas_data->nsampl = 0;
	meas_data->max_nsampl = max_nsampl;

	meas_data->sign = 0;

	meas_data->X_avg    = (double *)MKL_calloc(Norb, sizeof(double), MEM_DATA_ALIGN);
	meas_data->X_avg_sq = (double *)MKL_calloc(Norb, sizeof(double), MEM_DATA_ALIGN);
	meas_data->X_sq_avg = (double *)MKL_calloc(Norb, sizeof(double), MEM_DATA_ALIGN);
	meas_data->V_avg    = (double *)MKL_calloc(Norb, sizeof(double), MEM_DATA_ALIGN);
	meas_data->V_sq_avg = (double *)MKL_calloc(Norb, sizeof(double), MEM_DATA_ALIGN);
	meas_data->PE       = (double *)MKL_calloc(Norb, sizeof(double), MEM_DATA_ALIGN);
	meas_data->KE       = (double *)MKL_calloc(Norb, sizeof(double), MEM_DATA_ALIGN);
	meas_data->nX       = (double *)MKL_calloc(Norb, sizeof(double), MEM_DATA_ALIGN);

	meas_data->iteration_X_avg = (double *)MKL_calloc(Norb * max_nsampl, sizeof(double), MEM_DATA_ALIGN);
	meas_data->iteration_X0    = (double *)MKL_calloc(Norb * max_nsampl, sizeof(double), MEM_DATA_ALIGN);
	meas_data->n_local_accept = 0;
	meas_data->n_local_total = 0;
	meas_data->n_block_accept = 0;
	meas_data->n_block_total = 0;
	meas_data->n_flip_accept = 0;
	meas_data->n_flip_total = 0;
}


//________________________________________________________________________________________________________________________
///
/// \brief Delete phonon measurement data (free memory)
///
void DeletePhononData(measurement_data_phonon_t *restrict meas_data)
{
	MKL_free(meas_data->iteration_X0);
	MKL_free(meas_data->iteration_X_avg);
	MKL_free(meas_data->nX);
	MKL_free(meas_data->KE);
	MKL_free(meas_data->PE);
	MKL_free(meas_data->V_sq_avg);
	MKL_free(meas_data->V_avg);
	MKL_free(meas_data->X_sq_avg);
	MKL_free(meas_data->X_avg_sq);
	MKL_free(meas_data->X_avg);
}


//________________________________________________________________________________________________________________________
///
/// \brief Accumulate phonon measurements
///
void AccumulatePhononData(const greens_func_t *restrict Gu, const greens_func_t *restrict Gd, const int l_green, const double *restrict X, const double dt, const double *restrict omega, measurement_data_phonon_t *restrict meas_data)
{
	int l;
	int i;  // spatial index
	int o;  // orbital index

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

	double *cur_X_avg = (double *)MKL_calloc(Norb, sizeof(double), MEM_DATA_ALIGN);
	for (l = 0; l < L; l++)     // for all discrete time differences...
	{
		const int lplus = (l + 1) % L;
		for (o = 0; o < Norb; o++)
		{
			for (i = 0; i < Ncell; i++)
			{
				cur_X_avg[o] += signfac * X[i + o*Ncell + l*N];
				meas_data->X_sq_avg[o] += signfac * square(X[i + o*Ncell + l*N]);
				meas_data->iteration_X_avg[o + meas_data->nsampl * Norb] += signfac * X[i + o*Ncell + l*N];
				meas_data->V_avg[o] += signfac/dt * (X[i + o*Ncell + lplus*N] - X[i + o*Ncell + l*N]);
				meas_data->V_sq_avg[o] += signfac/(dt*dt) * square(X[i + o*Ncell + lplus*N] - X[i + o*Ncell + l*N]);
				meas_data->PE[o] += 0.5*square(omega[o])*signfac*square(X[i + o*Ncell + l*N]);
				meas_data->KE[o] -= 0.5/(dt*dt)*signfac*square(X[i + o*Ncell + lplus*N] - X[i + o*Ncell + l*N]);
			}
			meas_data->iteration_X0[o + meas_data->nsampl * Norb] += (sign/L) * X[o*Ncell + l*N];
		}
	}

	for (o = 0; o < Norb; o++)
	{
		meas_data->X_avg[o] += cur_X_avg[o];
		meas_data->X_avg_sq[o] += square(cur_X_avg[o]);
		meas_data->KE[0] += sign*0.5/dt;
	}

	// compute sum_i n_i X_i
	for (o = 0; o < Norb; o++)
	{
		double nX_avr = 0;
		for (i = 0; i < Ncell; i++)
		{
			const int io = i + Ncell * o;
			const double nu = (1 - Gu->mat[io + N*io]);
			const double nd = (1 - Gd->mat[io + N*io]);

			// phonon field at time slice of current Green's function
			nX_avr += (nu + nd) * X[i + o*Ncell + l_green*N];
		}
		nX_avr *= (sign / Ncell);

		meas_data->nX[o] += nX_avr;
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
	// total number of orbitals
	const int Norb  = meas_data->Norb;

	// normalization factor
	const double nfac = 1.0 / meas_data->nsampl;

	// divide by nsampl to get <s A>
	cblas_dscal(Norb, nfac, meas_data->X_avg, 1);
	cblas_dscal(Norb, nfac, meas_data->X_avg_sq, 1);
	cblas_dscal(Norb, nfac, meas_data->X_sq_avg, 1);
	cblas_dscal(Norb, nfac, meas_data->V_avg, 1);
	cblas_dscal(Norb, nfac, meas_data->V_sq_avg, 1);
	cblas_dscal(Norb, nfac, meas_data->PE, 1);
	cblas_dscal(Norb, nfac, meas_data->KE, 1);
	cblas_dscal(Norb, nfac, meas_data->nX, 1);

	// calculate average sign
	meas_data->sign *= nfac;
}


//________________________________________________________________________________________________________________________
///
/// \brief Print phonon measurement summary
///
void PrintPhononData(const measurement_data_phonon_t *meas_data)
{
	const double sgn = meas_data->sign;

	duprintf("_______________________________________________________________________________\n");
	duprintf("Summary of phonon data\n\n");
	duprintf("                    average sign: %g\n", sgn);
	int o;
	for (o = 0; o < meas_data->Norb; o++)
	{
		duprintf("\nResults for orbital %d\n", o);
		duprintf("          local acceptance ratio: %g\n", (double)meas_data->n_local_accept/meas_data->n_local_total);
		duprintf("          block acceptance ratio: %g\n", (double)meas_data->n_block_accept/meas_data->n_block_total);
		duprintf("           flip acceptance ratio: %g\n", (double)meas_data->n_flip_accept/meas_data->n_flip_total);
		duprintf("                       average X: %g\n", meas_data->X_avg[o]/sgn);
		duprintf("                     average X^2: %g\n", meas_data->X_sq_avg[o]/sgn);
		duprintf("                     average V^2: %g\n", meas_data->V_sq_avg[o]/sgn);
		duprintf("               average phonon KE: %g\n", meas_data->KE[o]/sgn);
		duprintf("               average phonon PE: %g\n", meas_data->PE[o]/sgn);
		duprintf("           average sum_i n_i X_i: %g\n", meas_data->nX[o]/sgn);
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
	sprintf(path, "%s_phonon_sign.dat",            fnbase); WriteData(path, &meas_data->sign,           sizeof(double), 1, false);
	sprintf(path, "%s_phonon_nsampl.dat",          fnbase); WriteData(path, &meas_data->nsampl,         sizeof(int),    1, false);
	sprintf(path, "%s_phonon_X_avg.dat",           fnbase); WriteData(path, meas_data->X_avg,           sizeof(double), Norb, false);
	sprintf(path, "%s_phonon_X_avg_sq.dat",        fnbase); WriteData(path, meas_data->X_avg_sq,        sizeof(double), Norb, false);
	sprintf(path, "%s_phonon_X_sq_avg.dat",        fnbase); WriteData(path, meas_data->X_sq_avg,        sizeof(double), Norb, false);
	sprintf(path, "%s_phonon_V_sq_avg.dat",        fnbase); WriteData(path, meas_data->V_sq_avg,        sizeof(double), Norb, false);
	sprintf(path, "%s_phonon_KE.dat",              fnbase); WriteData(path, meas_data->KE,              sizeof(double), Norb, false);
	sprintf(path, "%s_phonon_PE.dat",              fnbase); WriteData(path, meas_data->PE,              sizeof(double), Norb, false);
	sprintf(path, "%s_phonon_nX.dat",              fnbase); WriteData(path, meas_data->nX,              sizeof(double), Norb, false);

	sprintf(path, "%s_phonon_iteration_X_avg.dat", fnbase); WriteData(path, meas_data->iteration_X_avg, sizeof(double), Norb * max_nsampl, false);
	sprintf(path, "%s_phonon_iteration_X0.dat",    fnbase); WriteData(path, meas_data->iteration_X0,    sizeof(double), Norb * max_nsampl, false);
	sprintf(path, "%s_phonon_n_local_accept.dat",  fnbase); WriteData(path, &meas_data->n_local_accept, sizeof(int), 1, false);
	sprintf(path, "%s_phonon_n_local_total.dat",   fnbase); WriteData(path, &meas_data->n_local_total,  sizeof(int), 1, false);
	sprintf(path, "%s_phonon_n_block_accept.dat",  fnbase); WriteData(path, &meas_data->n_block_accept, sizeof(int), 1, false);
	sprintf(path, "%s_phonon_n_block_total.dat",   fnbase); WriteData(path, &meas_data->n_block_total,  sizeof(int), 1, false);
}
