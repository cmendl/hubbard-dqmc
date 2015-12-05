#include "progress.h"
#include <stdint.h>
#include <dupio.h>
#include <util.h>
#include <signal.h>

//________________________________________________________________________________________________________________________
///
/// \brief Store information on current progress of the DQMC simulation
///
typedef struct
{
	enum {INIT, EQUILIB, MEASURE, OUTPUT} state;	//!< current state
	int nequil, nsampl;								//!< total number of equilibration and measurement sweeps
	int niter;										//!< number of iterations in current state and current run
	int *iteration;									//!< pointer to number of total iterations (i.e. including previous checkpoints')
	uint64_t t_start;								//!< starting time (in nanoseconds)
	uint64_t t_equilib_start;						//!< starting time of warmup iterations
	uint64_t t_measure_start;						//!< starting time of measurement iterations
	uint64_t t_last_iter;							//!< time last iteration finished
}
dqmc_progress_t;

static dqmc_progress_t p;

static void PrintProgress(int signum)
{
	uint64_t t_current = GetTicks();
	duprintf("_______________________________________________________________________________\n");
	duprintf("Total elapsed time in current run (seconds): %.2f\n", (double)(t_current - p.t_start) / 1000000000.);
	switch (p.state)
	{
	case INIT:
		duprintf("Simulation is initializing.\n");
		break;
	case EQUILIB:
		duprintf("Completed equilibration iterations: %d out of %d\n", *p.iteration, p.nequil);
		if (p.niter > 0)
		{
			double time_per_iter = (double)(p.t_last_iter - p.t_equilib_start) / p.niter;
			duprintf("Estimated remaining time (equilibration): %.2f\n", (time_per_iter * (p.nequil - *p.iteration) + p.t_last_iter - t_current) / 1000000000.);
			duprintf("Estimated total remaining time (ignoring measurement cost): %.2f\n", (time_per_iter * (p.nequil + p.nsampl - *p.iteration) + p.t_last_iter - t_current) / 1000000000.);
		}
		break;
	case MEASURE:
		duprintf("Equilibration completed.\n");
		duprintf("Completed measurement iterations: %d out of %d\n", *p.iteration - p.nequil, p.nsampl);
		if (p.niter > 0)
		{
			double time_per_iter = (double)(p.t_last_iter - p.t_measure_start) / p.niter;
			duprintf("Estimated remaining time: %.2f\n", (time_per_iter * (p.nequil + p.nsampl - *p.iteration) + p.t_last_iter - t_current) / 1000000000.);
		}
		break;
	case OUTPUT:
		duprintf("Simulation completed, outputting data.\n");
		break;
	}

	#ifdef _WIN32
	signal(SIGBREAK, PrintProgress);
	#endif
}

int InitProgressTracking(int *iteration, const int nequil, const int nsampl)
{
	p.state = INIT;
	p.nequil = nequil;
	p.nsampl = nsampl;
	p.niter = 0;
	p.iteration = iteration;
	p.t_start = GetTicks();
	p.t_equilib_start = p.t_start;
	p.t_measure_start = p.t_start;
	p.t_last_iter = p.t_start;

	#ifdef _WIN32
	return (signal(SIGBREAK, PrintProgress) == SIG_ERR) ? -1 : 0;
	#else
	return sigaction(SIGUSR1, &(struct sigaction){.sa_handler = PrintProgress}, NULL);
	#endif
}

void UpdateProgress(void)
{
	uint64_t t_current = GetTicks();
	switch (p.state)
	{
	case INIT:
		p.state = (*p.iteration < p.nequil) ? EQUILIB : MEASURE;
		p.t_equilib_start = t_current;
		p.t_measure_start = t_current;
		break;
	case EQUILIB:
		p.t_last_iter = t_current;
		p.niter++;
		// + 1 because UpdateProgress is called before the for loop increments *p.iteration
		if (*p.iteration + 1 == p.nequil)
		{
			p.state = MEASURE;
			p.niter = 0;
			p.t_measure_start = t_current;
		}
		break;
	case MEASURE:
		p.t_last_iter = t_current;
		p.niter++;
		if (*p.iteration + 1 == p.nequil + p.nsampl)
		{
			p.state = OUTPUT;
		}
		break;
	case OUTPUT:
		break;
	}
}
