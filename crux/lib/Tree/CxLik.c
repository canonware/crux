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
	piSum += Pi[i];
    }
    if (piSum != 1.0) {
	for (i = 0; i < n; i++) {
	    Pi[i] /= piSum;
	}
    }

    // Fill in the upper triangle of Q.
    for (i = 0; i < n; i++) {
	Q[i*n + i] = 0.0;
    }
    for (i = 0; i < n; i++) {
	for (j = i + 1; j < n; j++) {
	    elm =  R[i*n + j] * Pi[j];
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

    // Q and V are operated on as column-major matrices in this function, in
    // order to interface cleanly with Fortran conventions.
    CxLikQ(n, Q, R, Pi);

    // Get the optimal work size.
    worksize = -1;
    dsyev_("V", "L", &n, Q, &n, qEigVals, &workQuery, &worksize, &info);
    CxmAssert(info == 0);
    worksize = workQuery;
    {
	double work[worksize];

	dsyev_("V", "L", &n, Q, &n, qEigVals, work, &worksize, &info);
	if (info != 0) {
	    return true;
	}
	// Q now contains the orthonormal eigenvectors of Q.
    }

    // Compute Q's inverse eigenvectors, and store the result in V.
    memcpy(V, Q, nSq * sizeof(double));
    info = clapack_dgetrf(CblasColMajor, n, n, V, n, ipiv);
    if (info != 0) {
	return true;
    }
    info = clapack_dgetri(CblasColMajor, n, V, n, ipiv);
    if (info != 0) {
	return true;
    }

    for (int i = 0; i < n; i++) {
	for (int j = 0; j < n; j++) {
	    for (int k = 0; k < n; k++) {
		qEigVecCube[i*nSq + j*n + k] = Q[k*n + i] * V[j*n + k];
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
	    for (int k = 0; k < n; k++) {
		p += qEigVecCube[i*nSq + j*n + k] * qEigValsExp[k];
	    }
	    if (p < 0.0) {
		p = 0.0;
	    }
	    P[i*n + j] = p;
	}
    }
}

// XXX Add support for Gamma-distributed rates.
CxmpInline void
CxLikExecuteStripe(CxtLik *lik, unsigned stripe) {
    double cL, stripeLnL;
    unsigned dim, cMin, cLim;
    double P[lik->dim * lik->dim];

    dim = lik->dim;
    cMin = lik->stripeWidth * stripe;
    cLim = cMin + lik->stripeWidth;

    // Iteratively process the execution plan.
    for (unsigned s = 0; s < lik->stepsLen; s++) {
	CxtLikStep *step = &lik->steps[s];
	CxLikPt(dim, P, step->model->qEigVecCube, step->model->qEigVals,
	  step->edgeLen);
	switch (step->variant) {
	    case CxeLikStepComputeCL: {
		for (unsigned c = cMin; c < cLim; c++) {
		    for (unsigned iP = 0; iP < dim; iP++) {
			cL = 0.0;
			for (unsigned iC = 0; iC < dim; iC++) {
			    cL += P[iP*dim + iC] * step->childMat[c*dim + iC];
			}
			step->parentMat[c*dim + iP] = cL;
		    }
		}
		break;
	    } case CxeLikStepMergeCL: {
		for (unsigned c = cMin; c < cLim; c++) {
		    for (unsigned iP = 0; iP < dim; iP++) {
			cL = 0.0;
			for (unsigned iC = 0; iC < dim; iC++) {
			    cL += P[iP*dim + iC] * step->childMat[c*dim + iC];
			}
			step->parentMat[c*dim + iP] *= cL;
		    }
		}
		break;
	    } default: {
		CxmNotReached();
	    }
	}
    }

    // For each model in the mixture, aggregrate conditional likelihoods and
    // weight them according to frequency priors, in order to compute
    // site-specific lnL's.
    for (unsigned m = 0; m < lik->modelsLen; m++) {
	CxtLikModel *model = &lik->models[m];
	if (model->weight == 0.0 || model->entire) {
	    // Model is not part of the mixture, or its lnL is already computed.
	    continue;
	}
	stripeLnL = 0.0;
	for (unsigned c = cMin; c < cLim; c++) {
	    double L = 0.0;
	    for (unsigned i = 0; i < dim; i++) {
		L += model->piDiag[i] * model->cL->mat[c*dim + i];
	    }
	    stripeLnL += log(L) * (double)lik->charFreqs[c];
	}
	model->stripeLnL[stripe] = stripeLnL;
    }
}

// XXX Generalize to handle threads.
bool
CxLikExecute(CxtLik *lik, double *rLnL) {
    unsigned stripe;
    double lnL;

    if (lik->stepsLen > 0) {
	// Compute log-likelihoods for each model in the mixture, stripewise.
	for (stripe = 0; stripe < lik->nstripes; stripe++) {
	    CxLikExecuteStripe(lik, stripe);
	}
    }

    // Sum stripeLnL for each model, and compute the weighted log-likelihood.
    lnL = 0.0;
    for (unsigned m = 0; m < lik->modelsLen; m++) {
	CxtLikModel *model = &lik->models[m];
	if (model->weight == 0.0) {
	    // Model is not part of the mixture.
	    continue;
	}
	if (model->entire == false) {
	    double modelLnL = 0.0;
	    for (unsigned s = 0; s < lik->nstripes; s++) {
		modelLnL += model->stripeLnL[s];
	    }
	    model->lnL = modelLnL;
	}
	lnL += model->lnL * model->weight;
    }
    *rLnL = lnL;

    return false;
}
