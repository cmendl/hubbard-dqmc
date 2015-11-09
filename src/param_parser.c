#include "param_parser.h"
#include "dupio.h"
#include "hash_table.h"
#include <mkl.h>
#include <stdlib.h>
#include <string.h>


//________________________________________________________________________________________________________________________
///
/// \brief Set default simulation parameters
///
void SetDefaultParameters(sim_params_t *params)
{
	// lattice field dimensions
	params->Nx = 6;
	params->Ny = 6;

	// Coulomb coupling constant in the Hubbard hamiltonian
	params->U = 4.0;

	// set t' (next-nearest neighbor) hopping parameter to zero
	params->tp = 0;

	// set chemical potential to zero (half-filling)
	params->mu = 0;

	// imaginary-time step size
	params->dt = 1.0/8;

	// number of time steps
	params->L = 16;
	// largest number of B_l matrices multiplied together before performing a QR decomposition
	params->prodBlen = 4;
	// number of "time slice wraps" before recomputing the Green's function
	params->nwraps = 8;

	// not taking phonons into account by default
	params->use_phonons = false;

	// fill phonon parameters even if 'use_phonons' is false by default
	params->phonon_params.omega = 1;
	params->phonon_params.g = 1;
	params->phonon_params.box_width = 8;
	params->phonon_params.nblock_updates = 0;	// disable block updates

	// UNIX time
	params->itime = time(NULL);
	//// artificial UNIX time
	//params->itime = 1420000241;

	// number of equilibration and sampling iterations
	params->nequil = 256;
	params->nsampl = 1024;

	// disable unequal time measurements by default
	params->nuneqlt = 0;
}


//________________________________________________________________________________________________________________________
///
/// \brief Parse parameter text file containing simulation parameters; parameters not set in the file remain unchanged
///
int ParseParameterFile(const char *filename, sim_params_t *params)
{
	FILE *fd = fopen(filename, "r");
	if (fd == NULL) {
		duprintf("Cannot open file '%s'.\n", filename);
		return -1;
	}

	// load the input file into a hash table
	ht_t params_table;
	htInit(&params_table, 64); // 64 buckets

	// read file line by line
	char line[1024];
	while (fgets(line, sizeof(line), fd) != NULL)
	{
		// skip text after a comment character
		strtok(line, "#");

		// delimiter characters
		const char *delim = " =:\t\r\n";

		// parameter name
		char *name = strtok(line, delim);
		if (name == NULL) {
			continue;
		}

		// parameter value
		char *value = strtok(NULL, delim);
		if (value == NULL)
		{
			duprintf("Missing value for parameter '%s' in input file '%s'.\n", name, filename);
			continue;
		}

		// add an entry to the hash table
		char *value_copy = (char *)MKL_malloc(strlen(value) + 1, MEM_DATA_ALIGN);
		strcpy(value_copy, value);

		if (htInsert(&params_table, name, value_copy) != 0)
		{
			duprintf("Warning: Duplicate parameter '%s' in input file '%s'.\n", name, filename);
		}
	}

	// update params with values from hash table
	char *value;
	if ((value = htGet(&params_table, "Nx"))       != NULL) params->Nx       = atoi(value);
	if ((value = htGet(&params_table, "Ny"))       != NULL) params->Ny       = atoi(value);
	if ((value = htGet(&params_table, "U"))        != NULL) params->U        = atof(value);
	if ((value = htGet(&params_table, "tp"))       != NULL) params->tp       = atof(value);
	if ((value = htGet(&params_table, "mu"))       != NULL) params->mu       = atof(value);
	if ((value = htGet(&params_table, "dt"))       != NULL) params->dt       = atof(value);
	if ((value = htGet(&params_table, "L"))        != NULL) params->L        = atoi(value);
	if ((value = htGet(&params_table, "prodBlen")) != NULL) params->prodBlen = atoi(value);
	if ((value = htGet(&params_table, "nwraps"))   != NULL) params->nwraps   = atoi(value);
	if ((value = htGet(&params_table, "use_phonons")) != NULL)
	{
		if (strcmp(value, "true") == 0) {
			params->use_phonons = true;
		}
		else if (strcmp(value, "false") == 0) {
			params->use_phonons = false;
		}
		else {
			duprintf("Warning: unrecognized 'use_phonons = %s' in file '%s', should be either 'true' or 'false'.\n", value, filename);
		}
	}
	if ((value = htGet(&params_table, "phonon_omega"))          != NULL) params->phonon_params.omega          = atof(value);
	if ((value = htGet(&params_table, "phonon_g"))              != NULL) params->phonon_params.g              = atof(value);
	if ((value = htGet(&params_table, "phonon_box_width"))      != NULL) params->phonon_params.box_width      = atof(value);
	if ((value = htGet(&params_table, "phonon_nblock_updates")) != NULL) params->phonon_params.nblock_updates = atof(value);
	if ((value = htGet(&params_table, "itime"))   != NULL) params->itime   = atoi(value);
	if ((value = htGet(&params_table, "nequil"))  != NULL) params->nequil  = atoi(value);
	if ((value = htGet(&params_table, "nsampl"))  != NULL) params->nsampl  = atoi(value);
	if ((value = htGet(&params_table, "nuneqlt")) != NULL) params->nuneqlt = atoi(value);

	// deallocate everything in the hash table
	htFree(&params_table);
	return 0;
}


//________________________________________________________________________________________________________________________
///
/// \brief Perform admissibility checks of the simulation parameters
///
int ValidateSimulationParameters(const sim_params_t *params)
{
	if (params->Nx <= 0 || params->Ny <= 0)
	{
		duprintf("Invalid parameters Nx = %i, Ny = %i: lattice dimensions must be positive.\n", params->Nx, params->Ny);
		return -1;
	}

	if (params->L <= 0)
	{
		duprintf("Invalid parameter L = %i: 'L' must be positive.\n", params->L);
		return -1;
	}

	if (params->L % params->prodBlen != 0 || params->L / params->prodBlen <= 0)
	{
		duprintf("Invalid parameters L = %i, prodBlen = %i: 'L' must be divisible by 'prodBlen' and 'L/prodBlen' must be positive.\n", params->L, params->prodBlen);
		return -1;
	}

	if (params->nwraps % params->prodBlen != 0 || params->nwraps/params->prodBlen <= 0)
	{
		duprintf("Invalid parameters nwraps = %i, prodBlen = %i: 'nwraps' must be divisible by 'prodBlen' and 'nwraps/prodBlen' must be positive.\n", params->nwraps, params->prodBlen);
		return -1;
	}

	return 0;
}

//________________________________________________________________________________________________________________________
///
/// \brief Print simulation parameters
///
void PrintParameters(const sim_params_t *params)
{
	duprintf("Simulation parameters\n\n");
	duprintf("               lattice dimension: %i x %i\n", params->Nx, params->Ny);
	duprintf("                               U: %g\n", params->U);
	duprintf("                              tp: %g\n", params->tp);
	duprintf("                              mu: %g\n", params->mu);
	duprintf("                       time step: %g\n", params->dt);
	duprintf("                               L: %i\n", params->L);
	duprintf("                            beta: %g\n", params->L * params->dt);	// inverse temperature
	duprintf("                        prodBlen: %i\n", params->prodBlen);
	duprintf("                          nwraps: %i\n", params->nwraps);
	duprintf("                   using phonons: %i\n", params->use_phonons);
	if (params->use_phonons)
	{
		duprintf("                phonon frequency: %g\n", params->phonon_params.omega);
		duprintf("               electron-phonon g: %g\n", params->phonon_params.g);
		duprintf("         phonon update box width: %g\n", params->phonon_params.box_width);
		duprintf("  number of phonon block updates: %i\n", params->phonon_params.nblock_updates);
	}
	duprintf("                           itime: %lli\n", params->itime);
	duprintf("        equilibration iterations: %i\n", params->nequil);
	duprintf("          measurement iterations: %i\n", params->nsampl);
	if (params->nuneqlt > 0)
	{
		duprintf("       unequal time measurements: every %i iteration(s)\n", params->nuneqlt);
	}
	else
	{
		duprintf("       unequal time measurements: disabled\n");
	}
	duprintf("\n");
}
