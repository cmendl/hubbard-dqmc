#ifndef MEASUREMENT_H
#define MEASUREMENT_H

#include "greens_func.h"


//________________________________________________________________________________________________________________________
///
/// \brief Measurement data structure
///
typedef struct
{
	double density_u;				//!< spin-up density
	double density_d;				//!< spin-down density
	double doubleocc;				//!< double occupancy

	int *latt_sum_map;				//!< lattice site index of coordinate sum of two lattice sites; matrix of size N x N

	double *uu_corr;				//!< up-up density correlations:         <n_{i,up} n_{j,up}>
	double *dd_corr;				//!< down-down density correlations:     <n_{i,dn} n_{j,dn}>
	double *ud_corr;				//!< up-down density cross correlations: <n_{i,up} n_{j,dn} + n_{i,dn} n_{j,up}>

	double *zz_corr;				//!< z-z spin correlations: <(n_{i,up} - n_{i,dn})(n_{j,up} - n_{j,dn})>
	double *xx_corr;				//!< x-x spin correlations: <(x_{i,+1} + x_{i,-1})(x_{j,+1} + x_{j,-1})> with x_{i,+1} = c^{dagger}_{i,dn} c_{i,up} and x_{i,-1} = c^{dagger}_{i,up} c_{i,dn}

	double sign;					//!< accumulated Green's function signs (+-1)
	int nsampl;						//!< number of accumulated samples

	int N;							//!< total number of lattice sites
}
measurement_data_t;


void AllocateMeasurementData(const int Nx, const int Ny, measurement_data_t *restrict meas_data);

void DeleteMeasurementData(measurement_data_t *restrict meas_data);


void AccumulateMeasurement(const greens_func_t *restrict Gu, const greens_func_t *restrict Gd, measurement_data_t *restrict meas_data);


void NormalizeMeasurementData(measurement_data_t *meas_data);


//________________________________________________________________________________________________________________________
///
/// \brief Unequal time measurement data structure
///
typedef struct
{
	double *Gtau0_u;			//!< concatenated unequal time spin-up   Green's functions G_u(tau,   0) with tau = 0, 1, ..., L-1; array of size L*N x N
	double *G0tau_u;			//!< concatenated unequal time spin-up   Green's functions G_u(0,   tau) with tau = 0, 1, ..., L-1; array of size N x L*N
	double *Geqlt_u;			//!< concatenated   equal time spin-up   Green's functions G_u(tau, tau) with tau = 0, 1, ..., L-1; array of size N x L*N
	double *Gtau0_d;			//!< concatenated unequal time spin-down Green's functions G_d(tau,   0) with tau = 0, 1, ..., L-1; array of size L*N x N
	double *G0tau_d;			//!< concatenated unequal time spin-down Green's functions G_d(0,   tau) with tau = 0, 1, ..., L-1; array of size N x L*N
	double *Geqlt_d;			//!< concatenated   equal time spin-down Green's functions G_d(tau, tau) with tau = 0, 1, ..., L-1; array of size N x L*N

	int *latt_sum_map;			//!< lattice site index of coordinate sum of two lattice sites; matrix of size N x N

	double *nn_corr;			//!< density correlations;  matrix of size N x L
	double *zz_corr;			//!< z-z spin correlations; matrix of size N x L
	double *xx_corr;			//!< x-x spin correlations; matrix of size N x L

	double *Hu;					//!< temporary matrix of size L*N x L*N for spin-up
	double *Hd;					//!< temporary matrix of size L*N x L*N for spin-down

	double sign;				//!< accumulated (equal time) Green's function signs (+-1)
	int nsampl;					//!< number of accumulated samples

	int N;						//!< total number of lattice sites
	int L;						//!< total number of time steps
}
measurement_data_unequal_time_t;


int AllocateUnequalTimeMeasurementData(const int Nx, const int Ny, const int L, measurement_data_unequal_time_t *restrict meas_data);

void DeleteUnequalTimeMeasurementData(measurement_data_unequal_time_t *restrict meas_data);


void AccumulateUnequalTimeMeasurement(const double sign, const double *const *Bu, const double *const *Bd, measurement_data_unequal_time_t *restrict meas_data);


void NormalizeUnequalTimeMeasurementData(measurement_data_unequal_time_t *meas_data);



#endif
