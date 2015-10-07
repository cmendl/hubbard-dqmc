#ifndef RANDOM_H
#define RANDOM_H

#include <stdint.h>


//________________________________________________________________________________________________________________________
///
/// \brief Random number generator based on George Marsaglia's Xorshift RNG algorithm
/// and improved by Sebastiano Vigna
///
typedef struct
{
	uint64_t s[16];		//!< seed consists of 16*64 = 1024 bits
	int p;				//!< counter
}
randseed_t;

void Random_SeedInit(uint64_t n, randseed_t *seed);


//________________________________________________________________________________________________________________________
//

uint64_t Random_GetUint(randseed_t *seed);

uint64_t Random_GetBoundedUint(randseed_t *seed, const uint64_t bound);

void Random_Shuffle(randseed_t *seed, const int n, int *seq);

double Random_GetUniform(randseed_t *seed);

double Random_GetNormal(randseed_t *seed);



#endif
