#include "random.h"
#include <math.h>
#include <assert.h>


//________________________________________________________________________________________________________________________
///
/// \brief Initialize seeds by applying Sebastiano Vigna's xorshift64* to the number 'n'
///
void Random_SeedInit(uint64_t n, randseed_t *seed)
{
	// must be non-zero
	assert(n > 0);

	int i;
	for (i = 0; i < 16; i++)
	{
		n ^= n >> 12;   // a
		n ^= n << 25;   // b
		n ^= n >> 27;   // c

		// 2685821657736338717 = 72821711 * 36882155347, from Pierre L'Ecuyer's paper
		seed->s[i] = n * 2685821657736338717LL;
	}

	seed->p = 0;
}


//________________________________________________________________________________________________________________________
///
/// \brief Use George Marsaglia's Xorshift RNG algorithm
/// improved by Sebastiano Vigna's xorshift1024* using a multiplication of the output
/// to draw an unsigned integer; generator period is (2^1024 - 1)
///
/// References:
///   - George Marsaglia. Xorshift RNGs. Journal of Statistical Software 8, 1-6 (2003)
///   - Sebastiano Vigna. An experimental exploration of Marsaglia's xorshift generators, scrambled. (2014)
///   - Pierre L'Ecuyer. Tables of linear congruential generators of different sizes and good lattice structure. Mathematics of Computation 68, 249-260 (1999) 
///
uint64_t Random_GetUint(randseed_t *seed)
{
	uint64_t s0 = seed->s[seed->p];
	seed->p = (seed->p + 1) & 15;
	uint64_t s1 = seed->s[seed->p];
	s1 ^= s1 << 31;     // a
	s1 ^= s1 >> 11;     // b
	s0 ^= s0 >> 30;     // c
	seed->s[seed->p] = s0 ^ s1;

	// 1181783497276652981 = 769 * 13611541 * 112902689, from Pierre L'Ecuyer's paper
	return seed->s[seed->p] * 1181783497276652981LL;
}


//________________________________________________________________________________________________________________________
///
/// \brief Draw a random unsigned integer from the interval [0, bound) without modulo bias,
/// see http://funloop.org/post/2015-02-27-removing-modulo-bias-redux.html
///
uint64_t Random_GetBoundedUint(randseed_t *seed, const uint64_t bound)
{
	const uint64_t thresh = (-bound) % bound;   // (-bound) = UINT64_MAX - bound + 1

	for (;;)
	{
		uint64_t r = Random_GetUint(seed);
		if (r >= thresh) {
			return r % bound;
		}
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief Random shuffle of (0, 1, ..., n-1) using the Fisher-Yates shuffle algorithm
///
void Random_Shuffle(randseed_t *seed, const int n, int *seq)
{
	int i;
	for (i = 0; i < n; i++)
	{
		int j = (int)Random_GetBoundedUint(seed, i + 1);
		if (j != i)
		{
			seq[i] = seq[j];
		}
		seq[j] = i;
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief Draw a uniform random sample from the interval (0, 1]
///
double Random_GetUniform(randseed_t *seed)
{
	// 0 <= u < 2^64
	uint64_t u = Random_GetUint(seed);
	// factor is 1 / 2^64
	return (u + 1.0) * 5.4210108624275221700372640e-20;
}


//________________________________________________________________________________________________________________________
///
/// \brief Draw a normal (Gaussian) random sample with mean 0 and standard deviation 1
///
double Random_GetNormal(randseed_t *seed)
{
	#define M_2PI 6.2831853071795864769

	// use Box-Muller transform
	double u1 = Random_GetUniform(seed);
	double u2 = Random_GetUniform(seed);
	double r = sqrt(-2*log(u1));
	return r * sin(M_2PI * u2);
}
