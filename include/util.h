#ifndef UTIL_H
#define UTIL_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <math.h>


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
