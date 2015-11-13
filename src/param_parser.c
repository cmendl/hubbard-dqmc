#include "param_parser.h"
#include "dupio.h"
#include "hash_table.h"
#include <mkl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


//________________________________________________________________________________________________________________________
///
/// \brief Reads a list like "1.2,3.4,5.6;7.8,9.0" into a double array. Commas and semi-colons are treated equivalently.
///
static int ReadList(const int Norb, double *a, char *list)
{
	char *value = strtok(list, ",;");
	int i;
	for (i = 0; i < Norb; i++)
	{
		// if value == NULL inside the loop, that means the list had too few entries
		if (value == NULL)
		{
			duprintf("ReadList() failed: too few entries in list.\n");
			return -1;
		}
		a[i] = atof(value);
		value = strtok(NULL, ",;");
	}

	// if value is not NULL here, that means the list had too many entries
	if (value != NULL)
	{
		duprintf("ReadList() failed: too many entries in list.\n");
		return -2;
	}
	return 0;
}


//________________________________________________________________________________________________________________________
///
/// \brief Reads a list like "0,2,3.456;1,5,3.141" into the array t. For this example,
///        we'd have t[0 + 2 * Norb] = 3.456 and t[1 + 5 * Norb] = 3.141.
///
static int ReadHoppingList(const int Norb, double *t, char *list)
{
	unsigned a, b; // orbital indices
	double tt;
	char *p = strtok(list, ";");
	while (p != NULL)
	{
		if (sscanf(p, " %u , %u , %lf ", &a, &b, &tt) != 3)
		{
			return -1; // bad format
		}
		if (a >= Norb || b >= Norb)
		{
			return -2; // orbital number out of bounds
		}
		// update the t array with bond information
		t[a + b * Norb] = tt;
		p = strtok(NULL, ";");
	}
	return 0;
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
			duprintf("Warning: Missing value for parameter '%s' in input file '%s'.\n", name, filename);
			continue;
		}

		// add the (name, value) pair to the hash table.
		// if hash table already has an entry for the parameter, append a semicolon and the new value.
		// this is so that multiple lines of t_ab or so can be concatenated.
		char *previous_value, *new_value;
		if ((previous_value = (char *)htGet(&params_table, name)) == NULL)
		{
			new_value = (char *)MKL_malloc(strlen(value) + 1, MEM_DATA_ALIGN);
			strcpy(new_value, value);
		}
		else
		{
			new_value = (char *)MKL_malloc(strlen(previous_value) + strlen(value) + 2, MEM_DATA_ALIGN);
			sprintf(new_value, "%s;%s", previous_value, value);
		}
		htInsert(&params_table, name, new_value); // this also frees previous_value
	}
	fclose(fd);

	// update params with values from hash table
	char *value;
	int Norb;

	if ((value = htGet(&params_table, "Norb")) != NULL)
	{
		Norb = atoi(value);
		if (Norb <= 0)
		{
			duprintf("Parameter 'Norb' must be positive.");
			return -2;
		}
		params->Norb = Norb;
		AllocateParameters(params);
	}
	else
	{
		duprintf("Parameter 'Norb' not found in input file '%s'.", filename);
		return -2;
	}
	if ((value = htGet(&params_table, "Nx")) != NULL) params->Nx = atoi(value);
	if ((value = htGet(&params_table, "Ny")) != NULL) params->Ny = atoi(value);

	if ((value = htGet(&params_table, "t_aa")) != NULL)
	{
		if (ReadHoppingList(Norb, params->t.aa, value) != 0)
		{
			duprintf("Error reading list for parameter 't_aa'.\n");
			return -3;
		}
	}
	if ((value = htGet(&params_table, "t_ab")) != NULL)
	{
		if (ReadHoppingList(Norb, params->t.ab, value) != 0)
		{
			duprintf("Error reading list for parameter 't_ab'.\n");
			return -3;
		}
	}
	if ((value = htGet(&params_table, "t_ac")) != NULL)
	{
		if (ReadHoppingList(Norb, params->t.ac, value) != 0)
		{
			duprintf("Error reading list for parameter 't_ac'.\n");
			return -3;
		}
	}
	if ((value = htGet(&params_table, "t_ad")) != NULL)
	{
		if (ReadHoppingList(Norb, params->t.ad, value) != 0)
		{
			duprintf("Error reading list for parameter 't_ad'.\n");
			return -3;
		}
	}
	if ((value = htGet(&params_table, "t_bc")) != NULL)
	{
		if (ReadHoppingList(Norb, params->t.bc, value) != 0)
		{
			duprintf("Error reading list for parameter 't_bc'.\n");
			return -3;
		}
	}
	if ((value = htGet(&params_table, "U")) != NULL)
	{
		if (ReadList(Norb, params->U, value) != 0)
		{
			duprintf("Error reading list for parameter 'U'.\n");
			return -3;
		}
	}
	if ((value = htGet(&params_table, "eps")) != NULL)
	{
		if (ReadList(Norb, params->eps, value) != 0)
		{
			duprintf("Error reading list for parameter 'eps'.\n");
			return -3;
		}
	}
	if ((value = htGet(&params_table, "mu")) != NULL) params->mu = atof(value);

	if ((value = htGet(&params_table, "phonon_omega")) != NULL)
	{
		if (ReadList(Norb, params->phonon_params.omega, value) != 0)
		{
			duprintf("Error reading list for parameter 'phonon_omega'.\n");
			return -3;
		}
	}
	if ((value = htGet(&params_table, "phonon_g")) != NULL)
	{
		if (ReadList(Norb, params->phonon_params.g, value) != 0)
		{
			duprintf("Error reading list for parameter 'phonon_g'.\n");
			return -3;
		}
	}
	if ((value = htGet(&params_table, "phonon_box_width")) != NULL) params->phonon_params.box_width = atof(value);
	if ((value = htGet(&params_table, "phonon_nblock_updates")) != NULL) params->phonon_params.nblock_updates = atof(value);

	params->use_phonons = false;
	int o;
	for (o = 0; o < Norb; o++)
	{
		if (params->phonon_params.g[o] != 0)
		{
			params->use_phonons = true;
			break;
		}
	}

	if ((value = htGet(&params_table, "dt"))       != NULL) params->dt       = atof(value);
	if ((value = htGet(&params_table, "L"))        != NULL) params->L        = atoi(value);
	if ((value = htGet(&params_table, "prodBlen")) != NULL) params->prodBlen = atoi(value);
	if ((value = htGet(&params_table, "nwraps"))   != NULL) params->nwraps   = atoi(value);

	if ((value = htGet(&params_table, "nequil"))   != NULL) params->nequil  = atoi(value);
	if ((value = htGet(&params_table, "nsampl"))   != NULL) params->nsampl  = atoi(value);
	if ((value = htGet(&params_table, "nuneqlt"))  != NULL) params->nuneqlt = atoi(value);
	if ((value = htGet(&params_table, "itime"))    != NULL) params->itime   = atoi(value);

	// deallocate everything in the hash table
	htFree(&params_table);
	return 0;
}


//________________________________________________________________________________________________________________________
///
/// \brief Allocate memory for parameters struct. Assumes params->Norb is already initialized.
///
void AllocateParameters(sim_params_t *params)
{
	const int Norb = params->Norb;
	params->t.aa                = (double *)MKL_calloc(Norb * Norb, sizeof(double), MEM_DATA_ALIGN);
	params->t.ab                = (double *)MKL_calloc(Norb * Norb, sizeof(double), MEM_DATA_ALIGN);
	params->t.ac                = (double *)MKL_calloc(Norb * Norb, sizeof(double), MEM_DATA_ALIGN);
	params->t.ad                = (double *)MKL_calloc(Norb * Norb, sizeof(double), MEM_DATA_ALIGN);
	params->t.bc                = (double *)MKL_calloc(Norb * Norb, sizeof(double), MEM_DATA_ALIGN);
	params->U                   = (double *)MKL_calloc(Norb,        sizeof(double), MEM_DATA_ALIGN);
	params->eps                 = (double *)MKL_calloc(Norb,        sizeof(double), MEM_DATA_ALIGN);
	params->phonon_params.omega = (double *)MKL_calloc(Norb,        sizeof(double), MEM_DATA_ALIGN);
	params->phonon_params.g     = (double *)MKL_calloc(Norb,        sizeof(double), MEM_DATA_ALIGN);
}


//________________________________________________________________________________________________________________________
///
/// \brief Free memory for parameters struct.
///
void DeleteParameters(sim_params_t *params)
{
	MKL_free(params->t.aa);
	MKL_free(params->t.ab);
	MKL_free(params->t.ac);
	MKL_free(params->t.ad);
	MKL_free(params->t.bc);
	MKL_free(params->U);
	MKL_free(params->eps);
	MKL_free(params->phonon_params.omega);
	MKL_free(params->phonon_params.g);
}


//________________________________________________________________________________________________________________________
///
/// \brief Perform admissibility checks of the simulation parameters
///
int ValidateSimulationParameters(const sim_params_t *params)
{
	const int Norb = params->Norb;
	if (Norb <= 0)
	{
		duprintf("Invalid parameter Norb = %i: must have at least one orbital.\n", params->Norb);
		return -1;
	}

	if (params->Nx <= 0 || params->Ny <= 0)
	{
		duprintf("Invalid parameters Nx = %i, Ny = %i: lattice dimensions must be positive.\n", params->Nx, params->Ny);
		return -1;
	}

	int i;
	for (i = 0; i < Norb; i++)
	{
		// within a unit cell, bonds must be between different orbitals
		if (params->t.aa[i + i*Norb] != 0.0)
		{
			duprintf("Invalid parameter t_aa: cannot have hopping between same orbital.\n");
			return -1;
		}

		// check for redundancy in intra-cell bonds
		int j;
		for (j = i + 1; j < Norb; j++)
		{
			if (params->t.aa[i + j * Norb] != 0.0 && params->t.aa[j + i * Norb] != 0.0)
			{
				duprintf("Invalid parameter t_aa: redundant intra-cell bond definitions.\n");
				return -1;
			}
		}
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


static void WriteList(const int Norb, const double *a, char *list)
{
	list[0] = '\0';
	char buf[32];
	int i;
	for (i = 0; i < Norb; i++)
	{
		snprintf(buf, sizeof(buf), "%g,", a[i]);
		strcat(list, buf);
	}
	// remove last comma
	list[strlen(list) - 1] = '\0';
}


static void WriteHoppingList(const int Norb, const double *t, char *list)
{
	list[0] = '\0';
	double tt;
	char buf[32];
	int i, j;
	for (i = 0; i < Norb; i++)
	{
		for (j = 0; j < Norb; j++)
		{
			tt = t[i + j * Norb];
			if (tt != 0.0)
			{
				snprintf(buf, sizeof(buf), "%d,%d,%g;", i, j, tt);
				strcat(list, buf);
			}
		}
	}
	// remove last semicolon
	list[strlen(list) - 1] = '\0';
}

//________________________________________________________________________________________________________________________
///
/// \brief Print simulation parameters
///
void PrintParameters(const sim_params_t *params)
{
	const int Norb = params->Norb;
	char list_buffer[512];

	duprintf("Simulation parameters\n\n");
	duprintf("              number of orbitals: %i\n", Norb);
	duprintf("               lattice dimension: %i x %i\n", params->Nx, params->Ny);
	WriteList(Norb, params->U, list_buffer);
	duprintf("                               U: %s\n", list_buffer);
	WriteList(Norb, params->eps, list_buffer);
	duprintf("                             eps: %s\n", list_buffer);
	WriteHoppingList(Norb, params->t.aa, list_buffer);
	duprintf("                            t_aa: %s\n", list_buffer);
	WriteHoppingList(Norb, params->t.ab, list_buffer);
	duprintf("                            t_ab: %s\n", list_buffer);
	WriteHoppingList(Norb, params->t.ac, list_buffer);
	duprintf("                            t_ac: %s\n", list_buffer);
	WriteHoppingList(Norb, params->t.ad, list_buffer);
	duprintf("                            t_ad: %s\n", list_buffer);
	WriteHoppingList(Norb, params->t.bc, list_buffer);
	duprintf("                            t_bc: %s\n", list_buffer);
	duprintf("                              mu: %g\n", params->mu);
	duprintf("                   using phonons: %i\n", params->use_phonons);
	if (params->use_phonons)
	{
		WriteList(Norb, params->phonon_params.omega, list_buffer);
		duprintf("                phonon frequency: %s\n", list_buffer);
		WriteList(Norb, params->phonon_params.g, list_buffer);
		duprintf("               electron-phonon g: %s\n", list_buffer);
		duprintf("         phonon update box width: %g\n", params->phonon_params.box_width);
		duprintf("  number of phonon block updates: %i\n", params->phonon_params.nblock_updates);
	}
	duprintf("                       time step: %g\n", params->dt);
	duprintf("                               L: %i\n", params->L);
	duprintf("                            beta: %g\n", params->L * params->dt);	// inverse temperature
	duprintf("                        prodBlen: %i\n", params->prodBlen);
	duprintf("                          nwraps: %i\n", params->nwraps);

	duprintf("                    initial seed: %lli\n", params->itime);
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
