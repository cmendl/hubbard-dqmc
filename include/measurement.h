#ifndef MEASUREMENT_H
#define MEASUREMENT_H

#include "greens_func.h"


//________________________________________________________________________________________________________________________
///
/// \brief Measurement data structure
///
typedef struct
{
	double density_u[2];			//!< spin-up density (mean value and average squares)
	double density_d[2];			//!< spin-down density (mean value and average squares)
	double doubleocc[2];			//!< double occupancy (mean value and average squares)

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


void AccumulateEqualTimeMeasurement(const greens_func_t *restrict Gu, const greens_func_t *restrict Gd, measurement_data_t *meas_data);


void NormalizeMeasurementData(measurement_data_t *meas_data);


//________________________________________________________________________________________________________________________
///
/// \brief Standard deviation given the mean and average square
///
static inline double SampleStandardDeviation(const int n, const double mean, const double sqr)
{
	// try to avoid numerical cancellation errors
	const double sqr2 = sqrt(sqr);
	return sqrt(n*(sqr2 - mean)*(sqr2 + mean) / (double)(n-1));
}



#endif
