#define CxMat_c
#include "CxMat.h"

// C-compatible prototype for the Fortran-based DSYEV in LAPACK.
extern void
dsyev_(char *jobz, char *uplo, int *n, double *A, int *lda, double *w,
  double *work, int *lwork, int *info);

// Decompose Q into eigenvectors/eigenvalues, then precompute a cube of
// products conceptually defined as:
//
//   eigVecCube[i][j][k] = eigVecs[i][k] * eigVecsInv[k][j]
//
//                    ___j
//                   /
//                  /  eigVecsInv
//                 i       |
//                   ______|_______
//      /|          /:     |      /
//     / |         / :     v     /|
//    j  i        /  :          / |
//               /   :         /  |
//   eigVecs----/.>  :        /   |
//             /_____________/    |
//       ____j |     ........|....|
//      /|     |    .        |    /
//     / |     |   .         |   /
//    /  |     |  .          |  /
//   k   i     | .           | /
//             |_____________|/
//  eigVecCube
//
// The above diagram shows how the two matrices can be oriented, such that each
// cube cell is the product of the projected matrix cells.
void
CxMatQDecomp(int n, double *Q, double *eigVecCube, double *eigVals) {
    int nSq = n*n;
    double U[nSq], V[nSq];
    double workQuery;
    int worksize, info;
    int ipiv[n];

    // The input matrix is destroyed, so make a copy of Q.
    memcpy(U, Q, n * n * sizeof(double));

    // Get the optimal work size.
    worksize = -1;
    dsyev_("V", "L", &n, U, &n, eigVals, &workQuery, &worksize, &info);
    CxmAssert(info == 0);
    worksize = workQuery;
    {
	double work[worksize];

	// From C's perspective, the upper triangle Q is valid, which makes it
	// the lower triangle from Fortran's perspective.
	dsyev_("V", "L", &n, U, &n, eigVals, work, &worksize, &info);
	CxmAssert(info == 0);
	// U now contains the orthonormal eigenvectors of Q.
    }

    // Compute U's inverse, and store the result in V.
    memcpy(V, U, n * n * sizeof(double));
    info = clapack_dgetri(CblasRowMajor, n, V, n, ipiv);
    CxmAssert(info == 0);

    for (int i = 0; i < n; i++) {
	for (int j = 0; j < n; j++) {
	    for (int k = 0; k < n; k++) {
		eigVecCube[i*nSq + j*n + k] = U[i*n + k] * V[k*n + j];
	    }
	}
    }
}

// Compute substitution probabilities, based on the eigenvector/eigenvalue
// decomposition of a Q matrix, mu, and t (muT == mu * t).
//
//           Q*mu*t
//   P(t) = e
void
CxMatPt(int n, double *P, double *eigVecCube, double *eigVals, double muT) {
    int nSq = n*n;
    double eigValsExp[n];
    double p;

    for (int i = 0; i < n; i++) {
	eigValsExp[i] = exp(eigVals[i] * muT);
    }

    for (int i = 0; i < n; i++) {
	for (int j = 0; j < n; j++) {
	    p = 0.0;
	    for (int x = 0; x < n; x++) {
		p += eigVecCube[i*nSq + j*n + x] * eigValsExp[x];
	    }
	    P[i*n + j] = p;
	}
    }
}
