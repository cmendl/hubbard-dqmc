#include "checkpoint.h"
#include <dupio.h>
#include <util.h>
#include <signal.h>
#include <sys/stat.h>
#include <stdint.h>


int stopped;

static void StopSimulation(int signum)
{
	duprintf("Signal %d received.\n", signum);
	stopped = 1;
}

int StopOnSIGINT(void)
{
	return (signal(SIGINT, StopSimulation) == SIG_ERR) ? -1 : 0;
	//return sigaction(SIGINT, &(struct sigaction){.sa_handler = StopSimulation}, NULL);
}


int SearchCheckpoint(const char *fnbase)
{
	char path[1024];
	struct stat st;

	sprintf(path, "%s_iteration.dat", fnbase);
	if (stat(path, &st) != 0) return -1;

	sprintf(path, "%s_randseed.dat", fnbase);
	if (stat(path, &st) != 0) return -1;

	sprintf(path, "%s_hsfield.dat", fnbase);
	if (stat(path, &st) != 0) return -1;

	return 0;
}


void LoadCheckpoint(const char *restrict fnbase, int *restrict iteration, randseed_t *restrict seed, spin_field_t *restrict s, const int LxN)
{
	char path[1024];
	sprintf(path, "%s_iteration.dat", fnbase);
	ReadData(path, iteration, sizeof(int), 1);

	sprintf(path, "%s_randseed.dat", fnbase);
	ReadData(path, seed->s, sizeof(uint64_t), 16);
	seed->p = 0;

	sprintf(path, "%s_hsfield.dat", fnbase);
	ReadData(path, s, sizeof(spin_field_t), LxN);
}


void SaveCheckpoint(const char *restrict fnbase, const int *restrict iteration, const randseed_t *restrict seed, const spin_field_t *restrict s, const int LxN)
{
	char path[1024];
	sprintf(path, "%s_iteration.dat", fnbase);
	WriteData(path, iteration, sizeof(int), 1, false);

	// instead of storing seed->p, just save a permutation of seed->s where the 0th element is the p'th element of seed->s.
	sprintf(path, "%s_randseed.dat", fnbase);
	uint64_t a[16];
	int i;
	for (i = 0; i < 16; i++)
	{
		a[i] = seed->s[(i - seed->p + 16) % 16];
	}
	WriteData(path, a, sizeof(uint64_t), 16, false);

	sprintf(path, "%s_hsfield.dat", fnbase);
	WriteData(path, s, sizeof(spin_field_t), LxN, false);
}
