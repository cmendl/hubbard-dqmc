#ifndef GREENS_FUNC_H
#define GREENS_FUNC_H

#include "time_flow.h"


//________________________________________________________________________________________________________________________
///
/// \brief Green's function matrix and corresponding determinant;
/// explicitly keep track of determinant since numerical computation based on matrix entries is not accurate enough
///
typedef struct
{
	double *mat;		//!< matrix entries, array of size N x N
	double logdet;		//!< logarithm of absolute value of determinant
	int    sgndet;		//!< sign of determinant
	int N;				//!< dimension
}
greens_func_t;


void AllocateGreensFunction(const int N, greens_func_t *G);

void DeleteGreensFunction(greens_func_t *G);

void CopyGreensFunction(const greens_func_t *restrict src, greens_func_t *restrict dst);


void GreenConstruct(const time_step_matrices_t *restrict tsm, const int slice_shift, greens_func_t *restrict G);

void GreenShermanMorrisonUpdate(const double delta, const int N, const int i, double *restrict Gmat);

void GreenTimeSliceWrap(const int N, const double *restrict B, const double *restrict invB, double *restrict Gmat);


void ComputeUnequalTimeGreensFunction(const int N, const int L, const time_step_matrices_t *restrict tsm, double *restrict H, double *restrict Gtau0, double *restrict G0tau, double *restrict Geqlt);



#endif
