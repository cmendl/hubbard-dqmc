#ifndef STRATONOVICH_H
#define STRATONOVICH_H

#include <stdint.h>


/// \brief Spin field type; zero encodes spin-up and one encodes spin-down
typedef uint8_t spin_field_t;


//________________________________________________________________________________________________________________________
///
/// \brief Parameters related to the Hubbard-Stratonovich spin field
///
typedef struct
{
	double *expVu[2];		//!< exp(-lambda), exp( lambda) with lambda = acosh(exp(dt*U/2)) for each orbital, used for spin-up   Green's function
	double *expVd[2];		//!< exp( lambda), exp(-lambda) with lambda = acosh(exp(dt*U/2)) for each orbital, used for spin-down Green's function
	double *delta[2];		//!< exp(2*lambda)-1, exp(-2*lambda)-1: delta parameter after a spin flip
	int Norb;				//!< number of orbitals per unit cell
}
stratonovich_params_t;


void FillStratonovichParameters(const int Norb, const double *U, const double dt, stratonovich_params_t *params);

void DeleteStratonovichParameters(stratonovich_params_t *params);



#endif
