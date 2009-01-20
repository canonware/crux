#include "CxLik.h"

#include <clapack.h>
#include <math.h>

// C-compatible prototype for the Fortran-based DSYEV in LAPACK.
extern void
dsyev_(char *jobz, char *uplo, int *n, double *A, int *lda, double *w,
  double *work, int *lwork, int *info);

// Compute the upper triangle of Q=R*Pi, and fill in the diagonals such that
// each row (of the corresponding full matrix) sums to 0.
CxmpInline void
CxLikQ(int n, double *Q, double *R, double *Pi) {
    double fixedRate, piSum, elm;
    unsigned i, j;

    // Rescale R, if necessary.
    fixedRate = R[(n-1)*n - 1];
    if (fixedRate != 1.0)
	for (i = 0; i < n; i++) {
	    for (j = i + 1; j < n; j++) {
		R[i*n + j] /= fixedRate;
	    }
	}

    // Rescale Pi, if necessary.
    piSum = 0.0;
    for (i = 0; i < n; i++) {
	piSum += Pi[i*n + i];
    }
    if (piSum != 1.0) {
	for (i = 0; i < n; i++) {
	    Pi[i*n + i] /= piSum;
	}
    }

    // Fill in the upper triangle of Q.
    for (i = 0; i < n; i++) {
	Q[i*n + i] = 0.0;
    }
    for (i = 0; i < n; i++) {
	for (j = i + 1; j < n; j++) {
	    elm =  R[i*n + j] * Pi[j*n + j];
	    Q[i*n + j] = elm;
	    Q[i*n + i] -= elm;
	    Q[j*n + j] -= elm;
	}
    }
}

// Decompose Q into eigenvectors/eigenvalues, then precompute a cube of
// products conceptually defined as:
//
//   qEigVecCube[i][j][k] = qEigVecs[i][k] * qEigVecsInv[k][j]
//
//                    ___j
//                   /
//                  / qEigVecsInv
//                 i       |
//                   ______|_______
//      /|          /:     |      /
//     / |         / :     v     /|
//    j  i        /  :          / |
//               /   :         /  |
//  qEigVecs----/.>  :        /   |
//             /_____________/    |
//       ____j |     ........|....|
//      /|     |    .        |    /
//     / |     |   .         |   /
//    /  |     |  .          |  /
//   k   i     | .           | /
//             |_____________|/
//  qEigVecCube
//
// The above diagram shows how the two matrices can be oriented, such that each
// cube cell is the product of the projected matrix cells.
bool
CxLikQDecomp(int n, double *R, double *Pi, double *qEigVecCube,
  double *qEigVals) {
    int nSq = n*n;
    double Q[nSq], V[nSq];
    double workQuery;
    int worksize, info;
    int ipiv[n];

    CxLikQ(n, Q, R, Pi);

    // Get the optimal work size.
    worksize = -1;
    dsyev_("V", "L", &n, Q, &n, qEigVals, &workQuery, &worksize, &info);
    CxmAssert(info == 0);
    worksize = workQuery;
    {
	double work[worksize];

	// From C's perspective, the upper triangle of Q is valid, which makes
	// it the lower triangle from Fortran's perspective.
	dsyev_("V", "L", &n, Q, &n, qEigVals, work, &worksize, &info);
	if (info != 0) {
	    return true;
	}
	// Q now contains the orthonormal eigenvectors of Q.
    }

    // Compute Q's inverse, and store the result in V.
    memcpy(V, Q, nSq * sizeof(double));
    info = clapack_dgetri(CblasRowMajor, n, V, n, ipiv);
    if (info != 0) {
	return true;
    }

    for (int i = 0; i < n; i++) {
	for (int j = 0; j < n; j++) {
	    for (int k = 0; k < n; k++) {
		qEigVecCube[i*nSq + j*n + k] = Q[i*n + k] * V[k*n + j];
	    }
	}
    }

    return false;
}

// Compute substitution probabilities, based on the eigenvector/eigenvalue
// decomposition of a Q matrix, mu, and t (muT == mu * t).
//
//           Q*mu*t
//   P(t) = e
void
CxLikPt(int n, double *P, double *qEigVecCube, double *qEigVals, double muT) {
    int nSq = n*n;
    double qEigValsExp[n];
    double p;

    for (int i = 0; i < n; i++) {
	qEigValsExp[i] = exp(qEigVals[i] * muT);
    }

    for (int i = 0; i < n; i++) {
	for (int j = 0; j < n; j++) {
	    p = 0.0;
	    for (int x = 0; x < n; x++) {
		p += qEigVecCube[i*nSq + j*n + x] * qEigValsExp[x];
	    }
	    P[i*n + j] = p;
	}
    }
}
