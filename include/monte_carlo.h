#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include "kinetic.h"
#include "greens_func.h"
#include "random.h"
#include "phonon.h"
#include "measurement.h"


void DQMCIteration(const kinetic_t *restrict kinetic, const stratonovich_params_t *restrict stratonovich_params, const int nwraps,
	randseed_t *restrict seed, spin_field_t *restrict s, time_step_matrices_t *restrict tsm_u, time_step_matrices_t *restrict tsm_d, greens_func_t *restrict Gu, greens_func_t *restrict Gd);

void DQMCPhononIteration(const double dt, const kinetic_t *restrict kinetic, const stratonovich_params_t *restrict stratonovich_params, const phonon_params_t *restrict phonon_params,
	const int nwraps, randseed_t *restrict seed, spin_field_t *restrict s, double *restrict X, double *restrict expX,
	time_step_matrices_t *restrict tsm_u, time_step_matrices_t *restrict tsm_d, greens_func_t *restrict Gu, greens_func_t *restrict Gd);


void PhononBlockUpdates(const double dt, const kinetic_t *restrict kinetic, const stratonovich_params_t *restrict stratonovich_params,
	const phonon_params_t *restrict phonon_params, randseed_t *restrict seed, const spin_field_t *restrict s, double *restrict X, double *restrict expX,
	time_step_matrices_t *restrict tsm_u, time_step_matrices_t *restrict tsm_d, greens_func_t *restrict Gu, greens_func_t *restrict Gd);


void DQMCSimulation(const double U, const double dt, const int L, const kinetic_t *restrict kinetic, const int prodBlen, const int nwraps,
	const int nequil, const int nsampl, randseed_t *restrict seed, measurement_data_t *restrict meas_data);

void DQMCPhononSimulation(const double U, const double dt, const int L, const kinetic_t *restrict kinetic, const int prodBlen, const int nwraps,
	const phonon_params_t *restrict phonon_params, const int nequil, const int nsampl, randseed_t *restrict seed, measurement_data_t *restrict meas_data);



#endif
