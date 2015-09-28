#ifndef GREENS_FUNC_H
#define GREENS_FUNC_H

#include "time_flow.h"


void GreenConstruct(const time_step_matrices_t *restrict tsm, const int slice_shift, double *restrict G);

void GreenShermanMorrisonUpdate(const double delta, const int N, const int i, double *restrict G);

void GreenTimeSliceWrap(const int N, const double *restrict B, const double *restrict invB, double *restrict G);



#endif
