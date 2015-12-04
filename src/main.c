#include "monte_carlo.h"
#include "kinetic.h"
#include "param_parser.h"
#include "random.h"
#include "stratonovich.h"
#include "profiler.h"
#include "checkpoint.h"
#include "progress.h"
#include "util.h"
#include "dupio.h"
#include <mkl.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <string.h>

// for sleep function and creating directories
#ifdef _WIN32
#include <windows.h>
#include <direct.h>

static int makedir(const char *path)
{
	return _mkdir(path);
}

#else

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

static int makedir(const char *path)
{
	return mkdir(path, 0755);
}

#endif

#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
#include <crtdbg.h>
#endif


int main(int argc, char *argv[])
{
	int status;

	// enable run-time memory check for debug builds
	#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	#endif

	// if compiled with OpenMP, use only 2 threads max for efficiency
	omp_set_num_threads(2);

	// initialize profiling
	Profile_Start();

	// make sure there's an input file
	if (argc < 2)
	{
		duprintf("No input file specified!.\n");
		return -1;
	}

	// zero simulation parameters and set initial time for random seed.
	sim_params_t params = { 0 };
	params.itime = GetTicks();

	// read parameters from input file.
	duprintf("Reading simulation parameters from file '%s'...\n", argv[1]);
	status = ParseParameterFile(argv[1], &params); // this also allocates memory for the arrays in params
	if (status < 0)
	{
		duprintf("Error parsing parameter file, exiting...\n");
		return -1;
	}

	// perform admissibility checks of the simulation parameters
	status = ValidateSimulationParameters(&params);
	if (status < 0) {
		duprintf("Parameters could not be validated, exiting...\n");
		return -2;
	}

	// create output directory or read from 2nd command line argument
	char path[1024];
	if (argc < 3)
	{
		duprintf("No output directory specified.\n");
		duprintf("Creating default directory based on parameters and initial seed...\n");
		makedir("output");
		sprintf(path, "output/N%ix%i_beta%g_mu%g_sim_%llu", params.Nx, params.Ny, params.L * params.dt, params.mu, params.itime);
	}
	else
	{
		strcpy(path, argv[2]);
	}
	makedir(path);
	duprintf("Using output directory '%s'.\n\n", path);
	
	// base output file name
	char fnbase[1024];
	sprintf(fnbase, "%s/sim_%llu", path, params.itime);

	// open simulation log file for writing
	sprintf(path, "%s_simulation.log", fnbase);
	fd_log = fopen(path, "a");
	if (fd_log == NULL)
	{
		duprintf("Cannot open log file '%s', exiting...\n", path);
		return -3;
	}

	// allocate and initialize equal time measurement data structure
	measurement_data_t meas_data;
	AllocateMeasurementData(params.Norb, params.Nx, params.Ny, &meas_data);

	// allocate and initialize unequal time measurement data structure
	measurement_data_unequal_time_t meas_data_uneqlt;
	if (params.nuneqlt > 0)
	{
		status = AllocateUnequalTimeMeasurementData(params.Norb, params.Nx, params.Ny, params.L, &meas_data_uneqlt);
		if (status < 0) {
			duprintf("Could not allocate unequal time measurement data structure (probably out of memory), exiting...\n");
			return -4;
		}
	}

	// state variables for DQMC simulation. these are also the variables stored in a checkpoint.
	int iteration = 0;
	randseed_t seed;
	const int LxN = params.L * params.Norb * params.Nx * params.Ny;
	spin_field_t *s = (spin_field_t *)MKL_malloc(LxN * sizeof(spin_field_t), MEM_DATA_ALIGN);

	// check for previous data in output directory
	if (SearchCheckpoint(fnbase) != 0) // no previous data found, start from scratch
	{
		// random generator seed; multiplicative constant from Pierre L'Ecuyer's paper
		Random_SeedInit(1865811235122147685LL * params.itime, &seed);

		// random initial Hubbard-Stratonovich field
		int i;
		for (i = 0; i < params.L * params.Norb * params.Nx * params.Ny; i++)
		{
			s[i] = (Random_GetUniform(&seed) < 0.5 ? 0 : 1);
		}

		// print a header and the parameters
		duprintf("Hubbard model DQMC\n------------------\n");
		duprintf("_______________________________________________________________________________\n");
		PrintSimulationParameters(&params);

		duprintf("Starting DQMC simulation...\n");
	}
	else // found a previous checkpoint, so load the previous state
	{
		LoadCheckpoint(fnbase, &iteration, &seed, s, LxN);
		duprintf("Loaded checkpoint.\n");
		if (iteration == params.nequil + params.nsampl)
		{
			duprintf("All iterations already completed in previous checkpoint. Exiting...\n");
			return -5;
		}

		// load previous measurement data
		LoadMeasurementData(fnbase, &meas_data);
		if (params.nuneqlt > 0)
		{
			LoadUnequalTimeMeasurementData(fnbase, &meas_data_uneqlt);
		}

		duprintf("Resuming DQMC simulation at iteration %d.\n", iteration);
	}

	// enable progress tracking (progress of simulation is shown whenever a SIGUSR1 signal is received)
	InitProgressTracking(&iteration, params.nequil, params.nsampl);

	// start timer
	const clock_t t_start = clock();

	// perform simulation
	DQMCSimulation(&params, &meas_data, &meas_data_uneqlt, &iteration, &seed, s);

	// stop timer
	const clock_t t_end = clock();
	double cpu_time = (double)(t_end - t_start) / CLOCKS_PER_SEC;
	duprintf("%d iterations completed, CPU time: %g\n", iteration, cpu_time);

	// at end of simulation, normalize measurement data and show summary of results
	if (iteration == params.nequil + params.nsampl)
	{
		duprintf("All iterations completed.\n");
		NormalizeMeasurementData(&meas_data);
		if (params.nuneqlt > 0)
		{
			NormalizeUnequalTimeMeasurementData(&meas_data_uneqlt);
		}
		// show some simulation results
		SummarizeMeasurementData(&meas_data);
	}

	// save checkpoint for next run. even if simulation is finished, saving the HS field
	// is helpful if someone wants to extend the simulation further.
	duprintf("Saving checkpoint to disk...");
	SaveCheckpoint(fnbase, &iteration, &seed, s, LxN);
	duprintf(" done.\n");

	// save simulation results as binary data to disk
	duprintf("Saving simulation results to disk...");
	SaveMeasurementData(fnbase, &meas_data);
	if (params.nuneqlt > 0)
	{
		SaveUnequalTimeMeasurementData(fnbase, &meas_data_uneqlt);
	}
	duprintf(" done.\n");

	// clean up
	Profile_Stop();
	fclose(fd_log);
	MKL_free(s);
	if (params.nuneqlt > 0)
	{
		DeleteUnequalTimeMeasurementData(&meas_data_uneqlt);
	}
	DeleteMeasurementData(&meas_data);
	DeleteSimulationParameters(&params);

	return stopped; // 0 if simulation ran to completion, 1 if stopped by SIGINT
}
