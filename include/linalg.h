#ifndef LINALG_H
#define LINALG_H


//________________________________________________________________________________________________________________________
///
/// \brief Compute the sign of 'x', can be -1, 0 or 1
///
static inline int Signum(const double x)
{
    return (0.0 < x) - (x < 0.0);
}


double Determinant(const int n, const double *restrict A);

int DeterminantSign(const int n, const double *restrict A);


int MatrixExp(const int n, const double *restrict A, double *restrict ret);


void MatrixProductSequence(const int n, const int num, const double **restrict A, double *restrict ret);


int BlockCyclicQR(const int n, const int p, double *restrict H, double *restrict tau);

int BlockCyclicTriangularInverse(const int n, const int p, double *restrict R);

int BlockCyclicInverseRotation(const int n, const int p, double *restrict A, const double *restrict tau);

int BlockCyclicInverse(const int n, const int p, double *restrict H);



#endif
