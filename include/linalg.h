#ifndef LINALG_H
#define LINALG_H

#define LAPACK_DISABLE_NAN_CHECK 

// wrapper for including CBLAS and LAPACKE

#ifdef USE_MKL

#include <mkl_cblas.h>
#include <mkl_lapacke.h>

#else

// generic include
#include <cblas.h>
#include <lapacke.h>

#endif


void domatcopy(size_t rows, size_t cols, const double alpha, const double *A, size_t lda, double *B, size_t ldb);


//________________________________________________________________________________________________________________________
//


int MatrixExp(const int n, const double *restrict A, double *restrict ret);


void MatrixProductSequence(const int n, const int num, const double *const *restrict A, double *restrict ret);


int BlockCyclicQR(const int n, const int p, double *restrict H, double *restrict tau);

int BlockCyclicTriangularInverse(const int n, const int p, double *restrict R);

int BlockCyclicInverseRotation(const int n, const int p, double *restrict A, const double *restrict tau);

int BlockCyclicInverse(const int n, const int p, double *restrict H);



#endif
