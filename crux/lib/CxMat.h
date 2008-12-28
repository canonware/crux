#ifndef CxMat_h
#define CxMat_h

#include "Cx.h"

#include <clapack.h>

#ifndef CxmUseInlines
double
CxMatDdet(unsigned n, double *A);
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
#endif

#endif // CxMat_h
