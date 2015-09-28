#ifndef KINETIC_H
#define KINETIC_H


//________________________________________________________________________________________________________________________
///
/// \brief Matrix exponential of the kinetic energy operator
///
typedef struct
{
	double *expK;			//!< exp(-dt K)
	double *inv_expK;		//!< exp( dt K)
	int N;					//!< total number of lattice sites
}
kinetic_t;


void NearestNeighborKineticExponential(const int Nx, const int Ny, const double mu, const double dt, kinetic_t *restrict kinetic);


void DeleteKineticExponential(kinetic_t *restrict kinetic);



#endif
