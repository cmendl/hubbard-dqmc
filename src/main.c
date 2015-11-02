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
	duprintf("                              tp: %g\n", params.tp);
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
	if (params.nuneqlt > 0)
	{
		duprintf("       unequal time measurements: every %i iteration(s)\n", params.nuneqlt);
	}
	else
	{
		duprintf("       unequal time measurements: disabled\n");
	}
	duprintf("\n");

	// calculate matrix exponential of the kinetic nearest neighbor hopping matrix
	kinetic_t kinetic;
	SquareLatticeKineticExponential(params.Nx, params.Ny, params.tp, params.mu, params.dt, &kinetic);

	// random generator seed; multiplicative constant from Pierre L'Ecuyer's paper
	randseed_t seed;
	Random_SeedInit(1865811235122147685LL * (uint64_t)params.itime, &seed);

	// allocate and initialize equal time measurement data structure
	measurement_data_t meas_data;
	AllocateMeasurementData(params.Nx, params.Ny, &meas_data);

	// allocate and initialize unequal time measurement data structure
	measurement_data_unequal_time_t meas_data_uneqlt;
	if (params.nuneqlt > 0)
	{
		status = AllocateUnequalTimeMeasurementData(params.Nx, params.Ny, params.L, &meas_data_uneqlt);
		if (status < 0) {
			duprintf("Could not allocate unequal time measurement data structure (probably out of memory), exiting...\n");
			return -3;
		}
	}

	duprintf("\nStarting DQMC simulation...\n");

	// start timer
	const clock_t t_start = clock();

	// perform simulation
	if (!params.use_phonons)
	{
		DQMCSimulation(params.U, params.dt, params.L, &kinetic, params.prodBlen, params.nwraps, params.nequil, params.nsampl, params.nuneqlt, &seed, &meas_data, &meas_data_uneqlt);
	}
	else
	{
		DQMCPhononSimulation(params.U, params.dt, params.L, &kinetic, params.prodBlen, params.nwraps, &params.phonon_params, params.nequil, params.nsampl, params.nuneqlt, &seed, &meas_data, &meas_data_uneqlt);
	}

	// normalize measurement data
	NormalizeMeasurementData(&meas_data);
	NormalizeUnequalTimeMeasurementData(&meas_data_uneqlt);

	// stop timer
	const clock_t t_end = clock();
	double cpu_time = (double)(t_end - t_start) / CLOCKS_PER_SEC;
	duprintf("Finished simulation, CPU time: %g\n", cpu_time);

	// show some simulation results
	duprintf("_______________________________________________________________________________\n");
	duprintf("Summary of simulation results\n\n");
	duprintf("                    average sign: %g\n", meas_data.sign);
	duprintf("           average total density: %g\n", meas_data.density_u + meas_data.density_d);
	duprintf("         average spin-up density: %g\n", meas_data.density_u);
	duprintf("       average spin-down density: %g\n", meas_data.density_d);
	duprintf("        average double occupancy: %g\n", meas_data.doubleocc);
	duprintf("            average local moment: %g\n", meas_data.density_u + meas_data.density_d - 2.0*meas_data.doubleocc);

	// save simulation results as binary data to disk
	duprintf("\nSaving simulation results to disk...");
	const int N = params.Nx * params.Ny;
	sprintf(path, "%s_density_u.dat", fnbase);	WriteData(path, &meas_data.density_u, sizeof(double), 1, false);
	sprintf(path, "%s_density_d.dat", fnbase);	WriteData(path, &meas_data.density_d, sizeof(double), 1, false);
	sprintf(path, "%s_doubleocc.dat", fnbase);	WriteData(path, &meas_data.doubleocc, sizeof(double), 1, false);
	sprintf(path, "%s_uu_corr.dat",   fnbase);	WriteData(path,  meas_data.uu_corr,   sizeof(double), N, false);
	sprintf(path, "%s_dd_corr.dat",   fnbase);	WriteData(path,  meas_data.dd_corr,   sizeof(double), N, false);
	sprintf(path, "%s_ud_corr.dat",   fnbase);	WriteData(path,  meas_data.ud_corr,   sizeof(double), N, false);
	sprintf(path, "%s_zz_corr.dat",   fnbase);	WriteData(path,  meas_data.zz_corr,   sizeof(double), N, false);
	sprintf(path, "%s_xx_corr.dat",   fnbase);	WriteData(path,  meas_data.xx_corr,   sizeof(double), N, false);
	sprintf(path, "%s_sign.dat",      fnbase);	WriteData(path, &meas_data.sign,      sizeof(double), 1, false);
	if (params.nuneqlt > 0)
	{
		sprintf(path, "%s_uneqlt_Gtau0_u.dat", fnbase);	WriteData(path, meas_data_uneqlt.Gtau0_u, sizeof(double), params.L*N*N, false);
		sprintf(path, "%s_uneqlt_G0tau_u.dat", fnbase);	WriteData(path, meas_data_uneqlt.G0tau_u, sizeof(double), params.L*N*N, false);
		sprintf(path, "%s_uneqlt_Geqlt_u.dat", fnbase);	WriteData(path, meas_data_uneqlt.Geqlt_u, sizeof(double), params.L*N*N, false);
		sprintf(path, "%s_uneqlt_Gtau0_d.dat", fnbase);	WriteData(path, meas_data_uneqlt.Gtau0_d, sizeof(double), params.L*N*N, false);
		sprintf(path, "%s_uneqlt_G0tau_d.dat", fnbase);	WriteData(path, meas_data_uneqlt.G0tau_d, sizeof(double), params.L*N*N, false);
		sprintf(path, "%s_uneqlt_Geqlt_d.dat", fnbase);	WriteData(path, meas_data_uneqlt.Geqlt_d, sizeof(double), params.L*N*N, false);

		sprintf(path, "%s_uneqlt_nn_corr.dat", fnbase);	WriteData(path, meas_data_uneqlt.nn_corr, sizeof(double), params.L*N, false);
		sprintf(path, "%s_uneqlt_zz_corr.dat", fnbase);	WriteData(path, meas_data_uneqlt.zz_corr, sizeof(double), params.L*N, false);
		sprintf(path, "%s_uneqlt_xx_corr.dat", fnbase);	WriteData(path, meas_data_uneqlt.xx_corr, sizeof(double), params.L*N, false);

		sprintf(path, "%s_uneqlt_sign.dat",    fnbase);	WriteData(path, &meas_data_uneqlt.sign,   sizeof(double), 1, false);
		sprintf(path, "%s_uneqlt_nsampl.dat",  fnbase);	WriteData(path, &meas_data_uneqlt.nsampl, sizeof(int),    1, false);
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
	DeleteKineticExponential(&kinetic);

	return 0;
}
