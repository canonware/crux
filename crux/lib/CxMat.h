#ifndef CxMat_h
#define CxMat_h

#include "Cx.h"

#include <clapack.h>
#include <math.h>

void
CxMatQDecomp(int n, double *Q, double *eigVecCube, double *eigVals);
void
CxMatPt(int n, double *P, double *eigVecCube, double *eigVals, double muT);

#ifndef CxmUseInlines
double
CxMatDdet(unsigned n, double *A);
double
CxMatLogDet(unsigned n, double *A);
#endif

#if (defined(CxmUseInlines) || defined(CxMat_c))
// Compute the determinant of an n x n square matrix.
CxmInline double
CxMatDdet(unsigned n, double *A) {
    double ret;
    int ipiv[n]; // Cython doesn't support this.
    int info;
    unsigned i;

    info = clapack_dgetrf(CblasRowMajor, n, n, A, n, ipiv);
    if (info != 0) {
	return 0.0;
    }

    ret = 1.0;
    for (i = 0; i < n; i++) {
	if (ipiv[i] != i) {
	    ret = -ret;
	}
	ret *= A[i * n + i];
    }

    return ret;
}

// Compute the LogDet/paralinear distance for two sequences, using A as an
// unscaled pairwise frequency matrix.
CxmInline double
CxMatLogDet(unsigned n, double *A) {
    double rowSums[n], colSums[n];
    double elm, detA, detRows, detCols, det;

    memset(rowSums, 0, sizeof(rowSums));
    memset(colSums, 0, sizeof(colSums));
    for (unsigned i = 0; i < n; i++) {
	for (unsigned j = 0; j < n; j++) {
	    elm = A[i*n + j];
	    rowSums[i] += elm;
	    colSums[j] += elm;
	}
    }

    detRows = detCols = 1.0;
    for (unsigned i = 0; i < n; i++) {
	detRows *= rowSums[i];
	detCols *= colSums[i];
    }
    if (detRows == 0.0 || detCols == 0.0) {
	return NAN;
    }

    // Do this after summing the rows/columns, since CxMatDdet() destroys A.
    detA = CxMatDdet(n, A);
    if (detA == 0.0) {
	return NAN;
    }

    det = detA / (sqrt(detRows) * sqrt(detCols));
    if (det == 1.0) {
	return 0.0;
    }

    return -log(det);
}
#endif

#endif // CxMat_h
