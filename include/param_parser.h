#ifndef PARAM_PARSER_H
#define PARAM_PARSER_H

#include "phonon.h"
#include <time.h>
#include <stdbool.h>


//________________________________________________________________________________________________________________________
///
/// \brief Simulation parameters
///
typedef struct
{
	int Nx;							//!< lattice field x dimension
	int Ny;							//!< lattice field y dimension

	double U;						//!< Coulomb coupling constant in the Hubbard hamiltonian
	double tp;						//!< t' (next-nearest neighbor) hopping parameter
	double mu;						//!< chemical potential in the Hubbard hamiltonian

	double dt;						//!< imaginary-time step size

	int L;							//!< number of time steps
	int prodBlen;					//!< largest number of B_l matrices multiplied together before performing a QR decomposition
	int nwraps;						//!< number of "time slice wraps" before recomputing the Green's function

	bool use_phonons;				//!< whether phonons should be included (Hubbard-Holstein model)
	phonon_params_t phonon_params;	//!< phonon parameters (only accessed if 'use_phonons' is true)

	time_t itime;					//!< UNIX time, used as seed for the random number generator

	int nequil;						//!< number of equilibration iterations
	int nsampl;						//!< number of sampling iterations
}
sim_params_t;


void SetDefaultParameters(sim_params_t *params);

int ParseParameterFile(const char *filename, sim_params_t *params);

int ValidateSimulationParameters(const sim_params_t *params);



#endif
