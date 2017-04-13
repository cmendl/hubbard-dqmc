#ifndef TIME_FLOW_H
#define TIME_FLOW_H

#include "kinetic.h"
#include "stratonovich.h"


//________________________________________________________________________________________________________________________
///
/// \brief Precomputed B matrices and products of these matrices for the fast evaluation of the time flow map
///
typedef struct
{
	double **B;             //!< array of B matrices for all time steps
	double **invB;          //!< array of inverse B matrices for all time steps
	double **Bprod;         //!< array of precomputed products of subsequent B matrices: (B_{prodBlen-1} ... B_1 B_0, B_{2prodlen-1} ... B_{prodBlen+1} B_{prodBlen}, ... )
	int N;                  //!< dimension (total number of lattice sites)
	int L;                  //!< total number of time steps
	int prodBlen;           //!< number of B matrices multiplied together in each entry of 'Bprod'; L must be a multiple of 'prodBlen'
	int numBprod;           //!< length of the array 'Bprod'; L = numBprod * prodBlen
}
time_step_matrices_t;


void AllocateTimeStepMatrices(const int N, const int L, const int prodBlen, time_step_matrices_t *restrict tsm);

void DeleteTimeStepMatrices(time_step_matrices_t *restrict tsm);

void CopyTimeStepMatrices(time_step_matrices_t *restrict src, time_step_matrices_t *restrict dst);

void InitTimeStepMatrices(const kinetic_t *restrict kinetic, const double *const expV[2], const spin_field_t *restrict s, time_step_matrices_t *restrict tsm);
void InitPhononTimeStepMatrices(const kinetic_t *restrict kinetic, const double *const expV[2], const spin_field_t *restrict s, const double *restrict expX, time_step_matrices_t *restrict tsm);

void UpdateTimeStepMatrices(const kinetic_t *restrict kinetic, const double *const expV[2], const spin_field_t *restrict s, const int l, time_step_matrices_t *restrict tsm);
void UpdatePhononTimeStepMatrices(const kinetic_t *restrict kinetic, const double *const expV[2], const spin_field_t *restrict s, const double *restrict expX, const int l, time_step_matrices_t *restrict tsm);


void TimeFlowMap(const time_step_matrices_t *restrict tsm, const int slice_shift, double *restrict Q, double *restrict tau, double *restrict d, double *restrict T);



#endif
