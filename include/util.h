#ifndef UTIL_H
#define UTIL_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <math.h>


#define MEM_DATA_ALIGN 64


#ifdef USE_MKL

#include <mkl_service.h>

#define algn_malloc(size) MKL_malloc(size, MEM_DATA_ALIGN)

#define algn_free MKL_free

#define algn_calloc(num, size) MKL_calloc(num, size, MEM_DATA_ALIGN)

#else

#include <malloc.h>
#ifdef __GNUC__
#include <mm_malloc.h>
#endif
#include <string.h>

#define algn_malloc(size) _mm_malloc(size, MEM_DATA_ALIGN)

#define algn_free _mm_free

static inline void *algn_calloc(size_t num, size_t size)
{
	size_t blocksize = num * size;
	void *ptr = algn_malloc(blocksize);
	if (ptr != NULL)
	{
		memset(ptr, 0, blocksize);
	}
	return ptr;
}

#endif


//________________________________________________________________________________________________________________________
///
/// \brief Square function x -> x^2
///
static inline double square(const double x)
{
	return x*x;
}


//________________________________________________________________________________________________________________________
///
/// \brief Uniform distance (infinity norm) between 'x' and 'y'
///
static inline double UniformDistance(const int n, const double *restrict x, const double *restrict y)
{
	double d = 0;
	int i;
	for (i = 0; i < n; i++)
	{
		d = fmax(d, fabs(x[i] - y[i]));
	}

	return d;
}


//________________________________________________________________________________________________________________________
//


int ReadData(const char *filename, void *data, const size_t size, const size_t n);

int WriteData(const char *filename, const void *data, const size_t size, const size_t n, const bool append);


//________________________________________________________________________________________________________________________
//


uint64_t GetTicks();

uint64_t GetTickRes();



#endif
