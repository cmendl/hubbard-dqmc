#ifndef LINALG_H
#define LINALG_H


int MatrixExp(const int n, const double *restrict A, double *restrict ret);


void MatrixProductSequence(const int n, const int num, const double *const *restrict A, double *restrict ret);


int BlockCyclicQR(const int n, const int p, double *restrict H, double *restrict tau);

int BlockCyclicTriangularInverse(const int n, const int p, double *restrict R);

int BlockCyclicInverseRotation(const int n, const int p, double *restrict A, const double *restrict tau);

int BlockCyclicInverse(const int n, const int p, double *restrict H);



#endif
