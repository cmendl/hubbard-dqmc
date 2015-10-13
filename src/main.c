#include "monte_carlo.h"
#include "kinetic.h"
#include "param_parser.h"
#include "profiler.h"
#include "util.h"
#include "dupio.h"
#include <mkl.h>
#include <math.h>

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

	Profile_Start();

	// (default) simulation parameters
	sim_params_t params;
	SetDefaultParameters(&params);

	// read parameters from file, overwriting default parameters
	if (argc >= 2) {
		duprintf("Reading simulation parameters from file '%s'...\n", argv[1]);
		status = ParseParameterFile(argv[1], &params);
		if (status < 0) {
			duprintf("Error parsing parameter file, exiting...\n");
			return -1;
		}
	}

	// perform admissibility checks of the simulation parameters
	status = ValidateSimulationParameters(&params);
	if (status < 0) {
		duprintf("Parameters could not be validated, exiting...\n");
		return -2;
	}

	// inverse temperature
	const double beta = params.L * params.dt;

	// trying to create output directory corresponding to parameters
	char path[1024];
	makedir("../output");
	sprintf(path, "../output/N%ix%i",                     params.Nx, params.Ny);                               makedir(path);
	sprintf(path, "../output/N%ix%i/U%g",                 params.Nx, params.Ny, params.U);                     makedir(path);
	sprintf(path, "../output/N%ix%i/U%g/beta%g",          params.Nx, params.Ny, params.U, beta);               makedir(path);
	sprintf(path, "../output/N%ix%i/U%g/beta%g/sim_%lli", params.Nx, params.Ny, params.U, beta, params.itime);
	while (makedir(path) < 0)
	{
		duprintf("Cannot create output directory '%s', changing 'itime'...\n", path);
		#ifdef _WIN32
		Sleep(15000);
		#else
		sleep(15);
		#endif
		params.itime += (clock() % 16) + 15;
		sprintf(path, "../output/N%ix%i/U%g/beta%g/sim_%lli", params.Nx, params.Ny, params.U, beta, params.itime);
	}
	duprintf("Created output directory '%s'...\n\n", path);

	// base output file name
	char fnbase[1024];
	sprintf(fnbase, "../output/N%ix%i/U%g/beta%g/sim_%lli/N%ix%i_U%g_beta%g_sim_%lli_ns%i", params.Nx, params.Ny, params.U, beta, params.itime, params.Nx, params.Ny, params.U, beta, params.itime, params.nsampl);

	// open simulation log file for writing
	sprintf(path, "%s_simulation.log", fnbase);
	fd_log = fopen(path, "w");

	duprintf("Hubbard model DQMC\n------------------\n");
	duprintf("_______________________________________________________________________________\n");
	duprintf("Simulation parameters\n\n");
	duprintf("               lattice dimension: %i x %i\n", params.Nx, params.Ny);
	duprintf("                               U: %g\n", params.U);
	duprintf("                              mu: %g\n", params.mu);
	duprintf("                       time step: %g\n", params.dt);
	duprintf("                               L: %i\n", params.L);
	duprintf("                            beta: %g\n", params.L * params.dt);	// inverse temperature
	duprintf("                        prodBlen: %i\n", params.prodBlen);
	duprintf("                          nwraps: %i\n", params.nwraps);
	duprintf("                   using phonons: %i\n", params.use_phonons);
	if (params.use_phonons)
	{
		duprintf("                phonon frequency: %g\n", params.phonon_params.omega);
		duprintf("               electron-phonon g: %g\n", params.phonon_params.g);
		duprintf("         phonon update box width: %g\n", params.phonon_params.box_width);
		duprintf("  number of phonon block updates: %i\n", params.phonon_params.nblock_updates);
	}
	duprintf("                           itime: %lli\n", params.itime);
	duprintf("        equilibration iterations: %i\n", params.nequil);
	duprintf("          measurement iterations: %i\n", params.nsampl);
	duprintf("\n");

	// calculate matrix exponential of the kinetic nearest neighbor hopping matrix
	kinetic_t kinetic;
	NearestNeighborKineticExponential(params.Nx, params.Ny, params.mu, params.dt, &kinetic);

	// random generator seed; multiplicative constant from Pierre L'Ecuyer's paper
	randseed_t seed;
	Random_SeedInit(1865811235122147685LL * (uint64_t)params.itime, &seed);

	// allocate and initialize measurement data structure
	measurement_data_t meas_data;
	AllocateMeasurementData(params.Nx, params.Ny, &meas_data);

	duprintf("\nStarting DQMC simulation...\n");

	// start timer
	const clock_t t_start = clock();

	// perform simulation
	if (!params.use_phonons)
	{
		DQMCSimulation(params.U, params.dt, params.L, &kinetic, params.prodBlen, params.nwraps, params.nequil, params.nsampl, &seed, &meas_data);
	}
	else
	{
		DQMCPhononSimulation(params.U, params.dt, params.L, &kinetic, params.prodBlen, params.nwraps, &params.phonon_params, params.nequil, params.nsampl, &seed, &meas_data);
	}

	// normalize measurement data
	NormalizeMeasurementData(&meas_data);

	// stop timer
	const clock_t t_end = clock();
	double cpu_time = (double)(t_end - t_start) / CLOCKS_PER_SEC;
	duprintf("Finished simulation, CPU time: %g\n", cpu_time);

	// show some simulation results
	duprintf("_______________________________________________________________________________\n");
	duprintf("Summary of simulation results\n\n");
	duprintf("                    average sign: %g\n", meas_data.sign);
	duprintf("         average spin-up density: %g +- %g\n", meas_data.density_u[0], SampleStandardDeviation(params.nsampl, meas_data.density_u[0], meas_data.density_u[1]));
	duprintf("       average spin-down density: %g +- %g\n", meas_data.density_d[0], SampleStandardDeviation(params.nsampl, meas_data.density_d[0], meas_data.density_d[1]));
	duprintf("        average double occupancy: %g +- %g\n", meas_data.doubleocc[0], SampleStandardDeviation(params.nsampl, meas_data.doubleocc[0], meas_data.doubleocc[1]));
	duprintf("            average local moment: %g\n", meas_data.density_u[0] + meas_data.density_d[0] - 2.0*meas_data.doubleocc[0]);

	// save simulation results as binary data to disk
	duprintf("\nSaving simulation results to disk...");
	sprintf(path, "%s_uu_corr.dat", fnbase);	WriteData(path, meas_data.uu_corr, sizeof(double), meas_data.N, false);
	sprintf(path, "%s_dd_corr.dat", fnbase);	WriteData(path, meas_data.dd_corr, sizeof(double), meas_data.N, false);
	sprintf(path, "%s_ud_corr.dat", fnbase);	WriteData(path, meas_data.ud_corr, sizeof(double), meas_data.N, false);
	sprintf(path, "%s_zz_corr.dat", fnbase);	WriteData(path, meas_data.zz_corr, sizeof(double), meas_data.N, false);
	sprintf(path, "%s_xx_corr.dat", fnbase);	WriteData(path, meas_data.xx_corr, sizeof(double), meas_data.N, false);
	sprintf(path, "%s_sign.dat",    fnbase);	WriteData(path, &meas_data.sign,   sizeof(double), 1,           false);
	duprintf(" done.\n");

	// clean up
	Profile_Stop();
	fclose(fd_log);
	DeleteKineticExponential(&kinetic);
	DeleteMeasurementData(&meas_data);

	return 0;
}
