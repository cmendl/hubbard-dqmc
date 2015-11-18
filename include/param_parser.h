#ifndef PARAM_PARSER_H
#define PARAM_PARSER_H

#include <stdint.h>
#include <stdbool.h>


//________________________________________________________________________________________________________________________
///
/// \brief Arrays for the hopping parameters
///
/// Consider four unit cells a, b, c, d
///
///  c    d
///  |\  /
///  | \/
///  | /\
///  |/  \
///  a----b
///
/// A bond_hoppings_t struct describes the minimum bonds between these unit cells such that
/// repeating the bonds everywhere on a rectangular lattice with periodic
/// boundary conditions tessellates all the bonds in the lattice.
///
/// Each t.xx element in the struct is an (Norb * Norb) length array.
/// For instance, if t.ab[0 + 2 * Norb] == 3.456, that means between orbital 0 of
/// unit cell a and orbital 2 of unit cell b, the hopping parameter is 3.456.
///
typedef struct
{
	double *aa;			//!< unit cell internal hopping
	double *ab;			//!< hopping between unit cells 'a' and 'b'
	double *ac;			//!< hopping between unit cells 'a' and 'c'
	double *ad;			//!< hopping between unit cells 'a' and 'd'
	double *bc;			//!< hopping between unit cells 'b' and 'c'
}
bond_hoppings_t;


void AllocateBondHoppings(const int Norb, bond_hoppings_t *bonds);

void DeleteBondHoppings(bond_hoppings_t *bonds);


//________________________________________________________________________________________________________________________
///
/// \brief Phonon parameters for the Hubbard-Holstein hamiltonian
///
typedef struct
{
	double *omega;				//!< phonon frequency, per orbital
	double *g;					//!< electron-phonon interaction strength, per orbital
	double box_width;			//!< delta updates of the phonon field variables are drawn from a uniform box-probability distribution on the interval [-box_width/2, box_width/2]
	int nblock_updates;			//!< number of phonon block updates at random sites, for each full iteration over lattice sites and time slices; set to zero to disable block updates
}
phonon_params_t;


//________________________________________________________________________________________________________________________
///
/// \brief Simulation parameters
///
typedef struct
{
	// lattice definition
	int Norb;						//!< number of orbitals per unit cell
	int Nx;							//!< lattice field x dimension
	int Ny;							//!< lattice field y dimension

	// Hamiltonian parameters
	bond_hoppings_t t;				//!< hopping parameters
	double *U;						//!< on-site repulsion, per orbital
	double *eps;					//!< site energies, per orbital
	double mu;						//!< chemical potential

	bool use_phonons;				//!< whether phonons should be included (Hubbard-Holstein model)
	phonon_params_t phonon_params;	//!< phonon parameters (only accessed if 'use_phonons' is true)

	// time flow parameters
	double dt;						//!< imaginary-time step size
	int L;							//!< number of time steps
	int prodBlen;					//!< largest number of B_l matrices multiplied together before performing a QR decomposition
	int nwraps;						//!< number of "time slice wraps" before recomputing the Green's function

	int nequil;						//!< number of equilibration iterations
	int nsampl;						//!< number of sampling iterations
	int nuneqlt;					//!< number of iterations before performing an unequal time measurement; set to 0 to disable unequal time measurements

	uint64_t itime;					//!< initial UNIX time, used as seed for the random number generator
}
sim_params_t;


int ParseParameterFile(const char *filename, sim_params_t *params);

void AllocateSimulationParameters(const int Norb, sim_params_t *params);

void DeleteSimulationParameters(sim_params_t *params);

int ValidateSimulationParameters(const sim_params_t *params);

void PrintSimulationParameters(const sim_params_t *params);



#endif
