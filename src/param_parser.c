#include "sim_params.h"
#include "dupio.h"
#include "hash_table.h"
#include <mkl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>


//________________________________________________________________________________________________________________________
///
/// \brief Allocate memory for bond hoppings
///
void AllocateBondHoppings(const int Norb, bond_hoppings_t *bonds)
{
	bonds->aa = (double *)MKL_calloc(Norb * Norb, sizeof(double), MEM_DATA_ALIGN);
	bonds->ab = (double *)MKL_calloc(Norb * Norb, sizeof(double), MEM_DATA_ALIGN);
	bonds->ac = (double *)MKL_calloc(Norb * Norb, sizeof(double), MEM_DATA_ALIGN);
	bonds->ad = (double *)MKL_calloc(Norb * Norb, sizeof(double), MEM_DATA_ALIGN);
	bonds->bc = (double *)MKL_calloc(Norb * Norb, sizeof(double), MEM_DATA_ALIGN);
}


//________________________________________________________________________________________________________________________
///
/// \brief Free memory for bond hoppings
///
void DeleteBondHoppings(bond_hoppings_t *bonds)
{
	MKL_free(bonds->bc);
	MKL_free(bonds->ad);
	MKL_free(bonds->ac);
	MKL_free(bonds->ab);
	MKL_free(bonds->aa);
}


//________________________________________________________________________________________________________________________
///
/// \brief Temporary structure for storing parameter values
///
typedef struct
{
	char *str[1024];		//!< string values
	int num;				//!< number of values
}
value_list_t;


//________________________________________________________________________________________________________________________
///
/// \brief Append parameter values to a list
///
static void AppendValues(value_list_t *list, const value_list_t *ap)
{
	int i;
	for (i = 0; i < ap->num; i++)
	{
		assert(0 <= list->num && list->num < 1024);
		list->str[list->num] = (char *)MKL_malloc((strlen(ap->str[i]) + 1) * sizeof(char), MEM_DATA_ALIGN);
		strcpy(list->str[list->num], ap->str[i]);
		list->num++;
	}
}

//________________________________________________________________________________________________________________________
///
/// \brief Delete parameter value list
///
static void DeleteValueList(value_list_t *list)
{
	int i;
	for (i = 0; i < list->num; i++)
	{
		MKL_free(list->str[i]);
	}

	list->num = 0;
}


//________________________________________________________________________________________________________________________
///
/// \brief Delete parameter value list contents and free structure memory
///
static void FreeValueListMemory(void *ptr)
{
	value_list_t *list = (value_list_t *)ptr;
	DeleteValueList(list);
	MKL_free(list);
}


//________________________________________________________________________________________________________________________
///
/// \brief Read a string list like { "0", "2", "3.4", "1", "5", "2.1" } into the array t;
/// for this example, we'd have t[0 + 2 * Norb] = 3.4 and t[1 + 5 * Norb] = 2.1.
///
static int ReadHoppings(const int Norb, const value_list_t *list, double *t)
{
	// must be triplets of the form i j d with i and j integers and d a floating-point number
	if (list->num % 3 != 0)
	{
		return -1;
	}

	int i;
	for (i = 0; i < list->num; i += 3)
	{
		 // orbital indices
		const int a = atoi(list->str[i]    );
		const int b = atoi(list->str[i + 1]);
		if (a < 0 || a >= Norb || b < 0 || b >= Norb)
		{
			return -2;	// orbital index out of bounds
		}

		// update the t array with bond information
		t[a + Norb*b] = atof(list->str[i + 2]);
	}

	return 0;
}


//________________________________________________________________________________________________________________________
///
/// \brief Reads a string list like { "1.2", "3.4", "5.6" } into a double array
///
static int ReadList(const int Norb, const value_list_t *list, double *a)
{
	// check if number of values matches
	if (list->num != Norb)
	{
		return -1;
	}

	int i;
	for (i = 0; i < Norb; i++)
	{
		a[i] = atof(list->str[i]);
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
	hash_table_t hashtable;
	AllocateHashTable(&hashtable, 64); // 64 buckets

	// read file line by line
	char line[1024];
	while (fgets(line, sizeof(line), fd) != NULL)
	{
		// skip text after a comment character
		char *csharp = strchr(line, '#');
		if (csharp != NULL)
		{
			// replace '#' by '\0';
			*csharp = '\0';
		}

		// delimiter characters
		const char *delim = " =:,;\t\r\n";

		// parameter name
		char *name = strtok(line, delim);
		if (name == NULL) {
			continue;
		}

		// parameter value(s)
		value_list_t value;
		value.num = 0;
		while ((value.str[value.num] = strtok(NULL, delim)) != NULL)
		{
			value.num++;
		}
		if (value.num == 0)
		{
			duprintf("Warning: Missing value for parameter '%s' in input file '%s'.\n", name, filename);
			continue;
		}

		// add the (name, value) pair to the hash table;
		// if the hash table already has an entry for the parameter name, concatenate values
		value_list_t *val_prev = (value_list_t *)HashTableGet(&hashtable, name);
		if (val_prev == NULL)	// not in hash table yet
		{
			value_list_t *v = MKL_calloc(1, sizeof(value_list_t), MEM_DATA_ALIGN);
			AppendValues(v, &value);
			HashTableInsert(&hashtable, name, v);
		}
		else
		{
			AppendValues(val_prev, &value);
		}
	}

	fclose(fd);

	// update parameters with values from hash table
	int Norb;
	value_list_t *value;
	if ((value = HashTableGet(&hashtable, "Norb")) != NULL)
	{
		assert(value->num > 0);
		Norb = atoi(value->str[0]);
		if (Norb <= 0)
		{
			duprintf("Parameter 'Norb' must be positive.");
			return -2;
		}
		AllocateSimulationParameters(Norb, params);
	}
	else
	{
		duprintf("Parameter 'Norb' not found in input file '%s'.", filename);
		return -2;
	}

	if ((value = HashTableGet(&hashtable, "Nx")) != NULL) { params->Nx = atoi(value->str[0]); }
	if ((value = HashTableGet(&hashtable, "Ny")) != NULL) { params->Ny = atoi(value->str[0]); }

	if ((value = HashTableGet(&hashtable, "pbc_shift")) != NULL) { params->pbc_shift = atoi(value->str[0]); }

	if ((value = HashTableGet(&hashtable, "t_aa")) != NULL)
	{
		if (ReadHoppings(Norb, value, params->t.aa) != 0)
		{
			duprintf("Error reading list for parameter 't_aa'.\n");
			return -3;
		}
	}
	if ((value = HashTableGet(&hashtable, "t_ab")) != NULL)
	{
		if (ReadHoppings(Norb, value, params->t.ab) != 0)
		{
			duprintf("Error reading list for parameter 't_ab'.\n");
			return -3;
		}
	}
	if ((value = HashTableGet(&hashtable, "t_ac")) != NULL)
	{
		if (ReadHoppings(Norb, value, params->t.ac) != 0)
		{
			duprintf("Error reading list for parameter 't_ac'.\n");
			return -3;
		}
	}
	if ((value = HashTableGet(&hashtable, "t_ad")) != NULL)
	{
		if (ReadHoppings(Norb, value, params->t.ad) != 0)
		{
			duprintf("Error reading list for parameter 't_ad'.\n");
			return -3;
		}
	}
	if ((value = HashTableGet(&hashtable, "t_bc")) != NULL)
	{
		if (ReadHoppings(Norb, value, params->t.bc) != 0)
		{
			duprintf("Error reading list for parameter 't_bc'.\n");
			return -3;
		}
	}
	if ((value = HashTableGet(&hashtable, "U")) != NULL)
	{
		if (ReadList(Norb, value, params->U) != 0)
		{
			duprintf("Error reading list for parameter 'U'.\n");
			return -3;
		}
	}
	if ((value = HashTableGet(&hashtable, "eps")) != NULL)
	{
		if (ReadList(Norb, value, params->eps) != 0)
		{
			duprintf("Error reading list for parameter 'eps'.\n");
			return -3;
		}
	}
	if ((value = HashTableGet(&hashtable, "mu")) != NULL) { params->mu = atof(value->str[0]); }

	// phonons

	if ((value = HashTableGet(&hashtable, "use_phonons")) != NULL)
	{
		assert(value->num > 0);

		if (strcmp(value->str[0], "true") == 0 || strcmp(value->str[0], "1") == 0) {
			params->use_phonons = true;
		}
		else if (strcmp(value->str[0], "false") == 0 || strcmp(value->str[0], "0") == 0) {
			params->use_phonons = false;
		}
		else {
			duprintf("Warning: unrecognized 'use_phonons = %s' in file '%s', should be either 'true'/'1' or 'false'/'0'.\n", value->str[0], filename);
		}
	}

	if ((value = HashTableGet(&hashtable, "phonon_omega")) != NULL)
	{
		if (ReadList(Norb, value, params->phonon_params.omega) != 0)
		{
			duprintf("Error reading list for parameter 'phonon_omega'.\n");
			return -3;
		}
	}
	if ((value = HashTableGet(&hashtable, "phonon_g")) != NULL)
	{
		if (ReadList(Norb, value, params->phonon_params.g) != 0)
		{
			duprintf("Error reading list for parameter 'phonon_g'.\n");
			return -3;
		}
	}
	if ((value = HashTableGet(&hashtable, "phonon_box_width")) != NULL)      { params->phonon_params.box_width      = atof(value->str[0]); }
	if ((value = HashTableGet(&hashtable, "phonon_nblock_updates")) != NULL) { params->phonon_params.nblock_updates = atoi(value->str[0]); }

	// time flow parameters
	if ((value = HashTableGet(&hashtable, "dt"))       != NULL) { params->dt       = atof(value->str[0]); }
	if ((value = HashTableGet(&hashtable, "L"))        != NULL) { params->L        = atoi(value->str[0]); }
	if ((value = HashTableGet(&hashtable, "prodBlen")) != NULL) { params->prodBlen = atoi(value->str[0]); }
	if ((value = HashTableGet(&hashtable, "nwraps"))   != NULL) { params->nwraps   = atoi(value->str[0]); }

	if ((value = HashTableGet(&hashtable, "nequil"))   != NULL) { params->nequil   = atoi(value->str[0]); }
	if ((value = HashTableGet(&hashtable, "nsampl"))   != NULL) { params->nsampl   = atoi(value->str[0]); }
	if ((value = HashTableGet(&hashtable, "neqlt"))    != NULL) { params->neqlt    = atoi(value->str[0]); }
	if ((value = HashTableGet(&hashtable, "nuneqlt"))  != NULL) { params->nuneqlt  = atoi(value->str[0]); }
	if ((value = HashTableGet(&hashtable, "itime"))    != NULL) { params->itime   = atoll(value->str[0]); }
	if ((value = HashTableGet(&hashtable, "max_time")) != NULL) { params->max_time = atoi(value->str[0]); }

	// deallocate everything in the hash table
	DeleteHashTable(&hashtable, FreeValueListMemory);

	return 0;
}


//________________________________________________________________________________________________________________________
///
/// \brief Allocate memory for the simulation parameters
///
void AllocateSimulationParameters(const int Norb, sim_params_t *params)
{
	params->Norb = Norb;

	AllocateBondHoppings(Norb, &params->t);

	params->U                   = (double *)MKL_calloc(Norb, sizeof(double), MEM_DATA_ALIGN);
	params->eps                 = (double *)MKL_calloc(Norb, sizeof(double), MEM_DATA_ALIGN);
	params->phonon_params.omega = (double *)MKL_calloc(Norb, sizeof(double), MEM_DATA_ALIGN);
	params->phonon_params.g     = (double *)MKL_calloc(Norb, sizeof(double), MEM_DATA_ALIGN);
}


//________________________________________________________________________________________________________________________
///
/// \brief Free memory for parameters struct
///
void DeleteSimulationParameters(sim_params_t *params)
{
	MKL_free(params->phonon_params.g);
	MKL_free(params->phonon_params.omega);
	MKL_free(params->eps);
	MKL_free(params->U);

	DeleteBondHoppings(&params->t);
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
		if (params->t.aa[i + i*Norb] != 0)
		{
			duprintf("Invalid parameter t_aa: cannot have hopping between same orbital %i.\n", i);
			return -1;
		}

		// check for redundancy in intra-cell bonds
		int j;
		for (j = i + 1; j < Norb; j++)
		{
			if (params->t.aa[i + j*Norb] != 0 && params->t.aa[j + i*Norb] != 0)
			{
				duprintf("Invalid parameter t_aa: redundant intra-cell %i <-> %i bond definitions.\n", i, j);
				return -1;
			}
		}
	}

	if (params->L <= 0)
	{
		duprintf("Invalid parameter L = %i: 'L' must be positive.\n", params->L);
		return -1;
	}

	if (params->prodBlen <= 0)
	{
		duprintf("Invalid parameter prodBlen = %i: 'prodBlen' must be positive.\n", params->prodBlen);
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

	if (params->itime == 0)
	{
		duprintf("Invalid parameter itime = 0; must be positive to be used as random seed.\n");
		return -1;
	}

	return 0;
}


//________________________________________________________________________________________________________________________
///
/// \brief Print the entries of a list
///
static inline void PrintList(const int num, const double *x)
{
	int i;
	for (i = 0; i < num-1; i++)
	{
		duprintf("%g, ", x[i]);
	}
	// avoid trailing comma
	duprintf("%g", x[num-1]);
}


//________________________________________________________________________________________________________________________
///
/// \brief Print the non-zero entries of a matrix as index-value tuples
///
static inline void PrintMatrix(const int num, const double *mat)
{
	int i, j;

	// print in row-major order
	for (i = 0; i < num; i++)
	{
		for (j = 0; j < num; j++)
		{
			const double x = mat[i + j*num];
			if (x != 0)
			{
				duprintf("(%d,%d,%g) ", i, j, x);
			}
		}
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief Print simulation parameters
///
void PrintSimulationParameters(const sim_params_t *params)
{
	const int Norb = params->Norb;

	duprintf("Simulation parameters\n\n");
	duprintf("              number of orbitals: %i\n", Norb);
	duprintf("               lattice dimension: %i x %i\n", params->Nx, params->Ny);
	duprintf("                       PBC shift: %i\n", params->pbc_shift);
	duprintf("                               U: ");   PrintList(Norb, params->U);       duprintf("\n");
	duprintf("                             eps: ");   PrintList(Norb, params->eps);     duprintf("\n");
	duprintf("                            t_aa: ");   PrintMatrix(Norb, params->t.aa);  duprintf("\n");
	duprintf("                            t_ab: ");   PrintMatrix(Norb, params->t.ab);  duprintf("\n");
	duprintf("                            t_ac: ");   PrintMatrix(Norb, params->t.ac);  duprintf("\n");
	duprintf("                            t_ad: ");   PrintMatrix(Norb, params->t.ad);  duprintf("\n");
	duprintf("                            t_bc: ");   PrintMatrix(Norb, params->t.bc);  duprintf("\n");
	duprintf("                              mu: %g\n", params->mu);
	duprintf("                   using phonons: %i\n", params->use_phonons);
	if (params->use_phonons)
	{
		duprintf("                phonon frequency: ");   PrintList(Norb, params->phonon_params.omega);  duprintf("\n");
		duprintf("               electron-phonon g: ");   PrintList(Norb, params->phonon_params.g);      duprintf("\n");
		duprintf("         phonon update box width: %g\n", params->phonon_params.box_width);
		duprintf("  number of phonon block updates: %i\n", params->phonon_params.nblock_updates);
	}
	duprintf("                       time step: %g\n", params->dt);
	duprintf("                               L: %i\n", params->L);
	duprintf("                            beta: %g\n", params->L*params->dt);	// inverse temperature
	duprintf("                        prodBlen: %i\n", params->prodBlen);
	duprintf("                          nwraps: %i\n", params->nwraps);

	duprintf("                    initial seed: %lli\n", params->itime);
	duprintf("        equilibration iterations: %i\n", params->nequil);
	duprintf("          measurement iterations: %i\n", params->nsampl);
	if (params->neqlt > 0)
	{
		duprintf("         equal time measurements: every %i time steps\n", params->neqlt);
	}
	else
	{
		duprintf("         equal time measurements: disabled\n");
	}
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
