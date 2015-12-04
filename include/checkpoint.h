#ifndef CHECKPOINT_H
#define CHECKPOINT_H

#include "random.h"
#include "stratonovich.h"


extern int stopped;


int StopOnSIGINT(void);


int InitCheckpointing(int nequil, int nsampl);


int LoadCheckpoint(const char *restrict fnbase, int *restrict iteration, randseed_t *restrict seed, spin_field_t *restrict s, const int LxN);


int SaveCheckpoint(const char *restrict fnbase, const int *restrict iteration, const randseed_t *restrict seed, const spin_field_t *restrict s, const int LxN);



#endif
