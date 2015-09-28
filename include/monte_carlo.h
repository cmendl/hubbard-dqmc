#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include "kinetic.h"
#include "random.h"
#include "phonon.h"
#include "measurement.h"


void DQMCSimulation(const double U, const double dt, const int L, const kinetic_t *restrict kinetic, const int prodBlen, const int nwraps,
	const int nequil, const int nsampl, randseed_t *restrict seed, measurement_data_t *restrict meas_data);

void DQMCPhononSimulation(const double U, const double dt, const int L, const kinetic_t *restrict kinetic, const int prodBlen, const int nwraps,
	const phonon_params_t *restrict phonon_params, const int nequil, const int nsampl, randseed_t *restrict seed, measurement_data_t *restrict meas_data);



#endif
