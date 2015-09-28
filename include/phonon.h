#ifndef PHONON_H
#define PHONON_H


//________________________________________________________________________________________________________________________
///
/// \brief Phonon parameters for the Hubbard-Holstein hamiltonian
///
typedef struct
{
	double omega;			//!< phonon frequency
	double g;				//!< electron-phonon interaction strength
	double box_width;		//!< delta updates of the phonon field variables are drawn from a uniform box-probability distribution on the interval [-box_width/2, box_width/2]
	int nblock_updates;		//!< number of phonon block updates at random sites, for each full iteration over lattice sites and time slices; set to zero to disable block updates
}
phonon_params_t;



#endif
