#include "CxLik.h"
#include "../CxMq.h"

#include <clapack.h>
#include <math.h>

// Worker thread context.
typedef struct {
    unsigned id;
    pthread_t pthread;
} CxtLikWorkerCtx;

// Message structure, passed through message queues.
typedef struct {
    CxtLik *lik;
    unsigned stripe;
} CxtLikMsg;

// Thread initialization control variable.
static pthread_once_t CxpLikOnce = PTHREAD_ONCE_INIT;

// Set to non-zero if threading is enabled, and initialization has occurred.
static unsigned CxpLikNThreads = 0;

// Array of CxpLikNThreads extant thread contexts.
static CxtLikWorkerCtx *CxpLikThreads;

// Message queues, used to communicate with worker threads.  Jobs are enqueued
// in CxpLikTodoMq, and the workers indicate job completion by sending the same
// message structures back through CxpLikDoneMq.
static CxtMq CxpLikTodoMq;
static CxtMq CxpLikDoneMq;

// C-compatible prototype for the Fortran-based DSYEV in LAPACK.
extern void
dsyev_(char *jobz, char *uplo, int *n, double *A, int *lda, double *w,
  double *work, int *lwork, int *info);

// Compute the upper triangle of Q=R*Pi, and fill in the diagonals such that
// each row (of the corresponding full matrix) sums to 0.
CxmpInline void
CxLikQ(int n, double *Q, double *R, double *Pi) {
    double fixedRate, rSum, rScale, piSum, elm;
    unsigned i, j;

    // Rescale R, if necessary.
    fixedRate = R[(n-1)*n - 1];
    if (fixedRate != 1.0) {
	for (i = 0; i < n; i++) {
	    for (j = i + 1; j < n; j++) {
		R[i*n + j] /= fixedRate;
	    }
	}
    }

    // Compute the scaling factor for R that adjusts the mean instantaneous
    // mutation rate to 1.
    rSum = 0.0;
    for (i = 0; i < n; i++) {
	for (j = i + 1; j < n; j++) {
	    rSum += R[i*n + j];
	}
    }
    rScale = ((double)n*n) / (rSum * 2.0);

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
	    elm =  R[i*n + j] * rScale * Pi[j];
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
	    P[i*n + j] = p;
	}
    }

    // Eigen decomposition can result in slightly negative probabilities, due
    // to rounding error.  Clamp such aberrant values in symmetric pairs, such
    // that slightly positive probabilities are also rounded to 0.0.
    for (int i = 0; i < n; i++) {
	for (int j = 0; j < n; j++) {
	    if (P[i*n + j] <= 0.0) {
		P[i*n + j] = 0.0;
		P[j*n + i] = 0.0;
	    }
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

// Worker thread entry function.
void *
CxpLikWorker(void *arg) {
//    CxtLikWorkerCtx *ctx = (CxtLikWorkerCtx *)arg;
    CxtLikMsg *msg;

    // Iteratively get a job, perform it, then send a return message to
    // indicate completion status.
    while (CxMqGet(&CxpLikTodoMq, &msg) == false) {
	CxLikExecuteStripe(msg->lik, msg->stripe);
	CxMqPut(&CxpLikDoneMq, &msg);
    }

    return NULL;
}

// Perform worker thread pool cleanup.
static void
CxpLikAtexit(void) {
    void *result;

    CxMqGetStop(&CxpLikTodoMq);
    for (unsigned i = 0; i < CxpLikNThreads; i++) {
	pthread_join(CxpLikThreads[i].pthread, &result);
    }
    free(CxpLikThreads);
    CxpLikThreads = NULL;
}

// Initialize the worker thread pool.  Since no errors are propagated from this
// function, perform initialization in an order that allows correct function
// (though degraded performance) even if an error occurs.
static void
CxpLikThreaded(void) {
    CxpLikThreads = (CxtLikWorkerCtx *)malloc(CxNcpus *
      sizeof(CxtLikWorkerCtx));
    if (CxpLikThreads == NULL) {
	return;
    }

    atexit(CxpLikAtexit);

    // Properly set the minimum number of message slots, so that the CxMq code
    // will never suffer from allocation failures after CxMqNew().
    if (CxMqNew(&CxpLikTodoMq, sizeof(CxtLikMsg *), CxNcpus)) {
	return;
    }
    if (CxMqNew(&CxpLikDoneMq, sizeof(CxtLikMsg *), CxNcpus)) {
	return;
    }

    for (unsigned i = 0; i < CxNcpus; i++) {
	CxpLikThreads[i].id = i;
	int err = pthread_create(&CxpLikThreads[i].pthread, NULL, CxpLikWorker,
	  (void *)&CxpLikThreads[i]);
	if (err) {
	    return;
	}
	CxpLikNThreads++;
    }
}

void
CxLikExecute(CxtLik *lik, double *rLnL) {
    double lnL;

    if (lik->stepsLen > 0) {
	if (CxNcpus > 1) {
	    pthread_once(&CxpLikOnce, CxpLikThreaded);
	}

	// Compute log-likelihoods for each model in the mixture, stripewise.
	if (CxpLikNThreads > 0 && lik->nstripes > 1) {
	    CxtLikMsg msgs[CxNcpus];
	    CxtLikMsg *msg;
	    unsigned stripe, lim, ndone;

	    // Stuff as many jobs as possible into the message queue.
	    if (lik->nstripes < sizeof(msgs) / sizeof(CxtLikMsg)) {
		lim = lik->nstripes;
	    } else {
		lim = sizeof(msgs) / sizeof(CxtLikMsg);
	    }
	    for (stripe = 0; stripe < lim; stripe++) {
		msgs[stripe].lik = lik;
		msgs[stripe].stripe = stripe;
		CxMqPut(&CxpLikTodoMq, &msgs[stripe]);
	    }

	    ndone = 0;
	    // Receive status messages and re-use the msg data structures to
	    // send more work messages.
	    for (; stripe < lik->nstripes; stripe++) {
		CxMqGet(&CxpLikDoneMq, &msg);
		msg->stripe = stripe;
		CxMqPut(&CxpLikTodoMq, msg);
		ndone++;
	    }

	    // Receive remaining status messages.
	    while (ndone < lik->nstripes) {
		CxMqGet(&CxpLikDoneMq, &msg);
		ndone++;
	    }
	} else {
	    // No worker threads; do all computations in the main thread.
	    for (unsigned stripe = 0; stripe < lik->nstripes; stripe++) {
		CxLikExecuteStripe(lik, stripe);
	    }
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
}
