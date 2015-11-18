#include "monte_carlo.h"
#include "kinetic.h"
#include "param_parser.h"
#include "profiler.h"
#include "util.h"
#include "dupio.h"
#include <mkl.h>
#include <math.h>
#include <omp.h>
#include <time.h>

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
		while (makedir(path) < 0)
		{
			duprintf("Cannot create output directory '%s', changing 'itime'...\n", path);
			#ifdef _WIN32
			Sleep(2000);
			#else
			sleep(2);
			#endif
			params.itime += (clock() % 16) + 15;
			sprintf(path, "output/N%ix%i_beta%g_mu%g_sim_%llu", params.Nx, params.Ny, params.L * params.dt, params.mu, params.itime);
		}
		duprintf("Created output directory '%s'.\n\n", path);
	}
	else
	{
		strcpy(path, argv[2]);
		duprintf("Using output directory '%s'.\n\n", path);
	}
	
	// base output file name
	char fnbase[1024];
	sprintf(fnbase, "%s/sim_%llu", path, params.itime);

	// open simulation log file for writing
	sprintf(path, "%s_simulation.log", fnbase);
	fd_log = fopen(path, "w");
	if (fd_log == NULL)
	{
		duprintf("Cannot open log file '%s'.\nCheck if output directory exists.\nexiting...\n", path);
		return -3;
	}

	duprintf("Hubbard model DQMC\n------------------\n");
	duprintf("_______________________________________________________________________________\n");
	PrintParameters(&params);

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

	duprintf("\nStarting DQMC simulation...\n");

	// start timer
	const clock_t t_start = clock();

	// perform simulation
	DQMCSimulation(&params, &meas_data, &meas_data_uneqlt);

	// normalize measurement data
	NormalizeMeasurementData(&meas_data);
	if (params.nuneqlt > 0)
	{
		NormalizeUnequalTimeMeasurementData(&meas_data_uneqlt);
	}

	// stop timer
	const clock_t t_end = clock();
	double cpu_time = (double)(t_end - t_start) / CLOCKS_PER_SEC;
	duprintf("Finished simulation, CPU time: %g\n", cpu_time);

	// show some simulation results
	duprintf("_______________________________________________________________________________\n");
	duprintf("Summary of simulation results\n\n");
	duprintf("                    average sign: %g\n", meas_data.sign);
	double total_density = 0.0;
	int i;
	for (i = 0; i < params.Norb; i++)
	{
		total_density += meas_data.density_u[i] + meas_data.density_d[i];
	}
	duprintf("           average total density: %g\n", total_density);
	for (i = 0; i < params.Norb; i++)
	{
		duprintf("\nResults for orbital %d\n", i);
		duprintf("           average total density: %g\n", meas_data.density_u[i] + meas_data.density_d[i]);
		duprintf("         average spin-up density: %g\n", meas_data.density_u[i]);
		duprintf("       average spin-down density: %g\n", meas_data.density_d[i]);
		duprintf("        average double occupancy: %g\n", meas_data.doubleocc[i]);
		duprintf("            average local moment: %g\n", meas_data.density_u[i] + meas_data.density_d[i] - 2.0*meas_data.doubleocc[i]);
	}

	// save simulation results as binary data to disk
	duprintf("\nSaving simulation results to disk...");
	SaveMeasurementData(fnbase, &meas_data);
	if (params.nuneqlt > 0)
	{
		SaveUnequalTimeMeasurementData(fnbase, &meas_data_uneqlt);
	}
	duprintf(" done.\n");

	// clean up
	Profile_Stop();
	fclose(fd_log);
	if (params.nuneqlt > 0)
	{
		DeleteUnequalTimeMeasurementData(&meas_data_uneqlt);
	}
	DeleteMeasurementData(&meas_data);
	DeleteParameters(&params);

	return 0;
}
