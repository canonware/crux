#include "CxLik.h"
#include "../CxMq.h"
#include "../CxLapack.h"

#include <math.h>

//#define CxmLikDebug

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

CxmpInline unsigned
CxpLikNxy2i(unsigned n, unsigned x, unsigned y) {
    CxmAssert(x < n);
    CxmAssert(y < n);
    CxmAssert(x < y);

#if 1
    CxmAssert(x != y);
#else
    if (x > y) {
	unsigned t = x;
	x = y;
	y = t;
    }
#endif

    return n*x + y - (((x+3)*x) >> 1) - 1;
}

// Compute Q=R*Pi, and fill in the diagonals such that each row (of the
// corresponding full matrix) sums to 0.  Q is a full matrix, R is an upper
// triangle vector, and Pi is a diagonal vector.  Also compute qNorm, which
// is a normalization factor that adjusts the mean substitution rate to 1.
CxmpInline void
CxLikQ(int n, double *qNorm, double *Q, double *RTri, double *PiDiag,
  double *PiDiagNorm) {
    double piSum, elm, qScale;
    unsigned i, j;

#ifdef CxmLikDebug
    fprintf(stderr, "R:\n");
    for (unsigned i = 0; i < n; i++) {
	for (unsigned j = 0; j < n; j++) {
	    if (i == j) {
		fprintf(stderr, " %10s", "-");
	    } else if (i < j) {
		fprintf(stderr, " %10.6f", RTri[CxpLikNxy2i(n,i,j)]
		  / RTri[CxpLikNxy2i(n,n-2,n-1)]);
	    } else {
		fprintf(stderr, " %10.6f", RTri[CxpLikNxy2i(n,j,i)]
		  / RTri[CxpLikNxy2i(n,n-2,n-1)]);
	    }
	}
	fprintf(stderr, "\n");
    }
#endif

    // Compute the normalization factor for Pi that adjusts the sum of Pi to 1,
    // and store the normalized Pi for later use.
    piSum = 0.0;
    for (i = 0; i < n; i++) {
	piSum += PiDiag[i];
    }
    for (i = 0; i < n; i++) {
	PiDiagNorm[i] = PiDiag[i] / piSum;
    }
#ifdef CxmLikDebug
    fprintf(stderr, "Pi:\n");
    for (unsigned i = 0; i < n; i++) {
	for (unsigned j = 0; j < n; j++) {
	    if (i == j) {
		fprintf(stderr, " %10.6f", PiDiagNorm[i]);
	    } else {
		fprintf(stderr, " %10s", "-");
	    }
	}
	fprintf(stderr, "\n");
    }
#endif

    // Fill in Q, in column-major form, and compute qScale (inverse of qNorm).
    qScale = 0.0;
    for (i = 0; i < n; i++) {
	for (j = 0; j < i; j++) {
	    elm = RTri[CxpLikNxy2i(n,j,i)] * PiDiagNorm[j];
	    qScale += elm * PiDiagNorm[i] * 2.0;
	    Q[i + j*n] = elm;
	}
	for (j = i + 1; j < n; j++) {
	    elm = RTri[CxpLikNxy2i(n,i,j)] * PiDiagNorm[j];
	    Q[i + j*n] = elm;
	}
    }
    // Set the diagonal such that each row sums to 0.
    for (i = 0; i < n; i++) {
	double diag = 0.0;
	for (j = 0; j < i; j++) {
	    diag -= Q[i + j*n];
	}
	for (j = i + 1; j < n; j++) {
	    diag -= Q[i + j*n];
	}
	Q[i + i*n] = diag;
    }

    // Return the normalization factor for Q.
    *qNorm = 1.0 / qScale;
#ifdef CxmLikDebug
    fprintf(stderr, "qNorm: %10.6f\n", *qNorm);
    fprintf(stderr, "Q:\n");
    for (unsigned i = 0; i < n; i++) {
	for (unsigned j = 0; j < n; j++) {
	    fprintf(stderr, " %10.6f", Q[i + j*n]);
	}
	fprintf(stderr, "\n");
    }
#endif
}

// Decompose Q into eigenvectors/eigenvalues, then precompute a cube of
// products conceptually defined as:
//
//   qEigVecCube[i][j][k] = qEigVecs[i][k] * qEigVecsInv[k][j]
//
// This cube is used when exponentiating Q, based on the fact that
//
//    Q      D  -1
//   e  = V*e *V  ,
//
// where V is the matrix of eigenvectors, and D is the diagonal matrix of
// eigenvalues.
//
//                    ___j
//                   /
//                  / qEigVecsInv
//                 i       |
//                   ______|_______
//      /|          /:     |      /|
//     / |         / :     v     / |
//    j  i        /  :          /  |
//               /   :         /   |
//  qEigVecs----/.>  :        /    |
//             /_____________/     |
//       ____j |     ........|.....|
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
CxLikQDecomp(int n, double *RTri, double *PiDiag, double *PiDiagNorm,
  double *qNorm, double *qEigVecCube, double *qEigVals) {
    int nSq = n*n;
    double Q[nSq], QEigVecs[nSq], V[nSq];
    double qEigValsI[n];
    double workQuery;
    int worksize, info;
    int ipiv[n];

    // Q and V are operated on as column-major matrices in this function, in
    // order to interface cleanly with Fortran conventions.
    CxLikQ(n, qNorm, Q, RTri, PiDiag, PiDiagNorm);

    // Get the optimal work size.
    worksize = -1;

    dgeev_("N", "V", &n, Q, &n, qEigVals, qEigValsI, NULL, &n, QEigVecs, &n,
      &workQuery, &worksize, &info);
    CxmAssert(info == 0);
    worksize = workQuery;
    {
	double work[worksize];

	dgeev_("N", "V", &n, Q, &n, qEigVals, qEigValsI, NULL, &n, QEigVecs,
	  &n, work, &worksize, &info);
	if (info != 0) {
	    return true;
	}
    }

    // Compute Q's inverse eigenvectors, and store the result in V.
    memcpy(V, QEigVecs, nSq * sizeof(double));
    dgetrf_(&n, &n, V, &n, ipiv, &info);
    if (info != 0) {
	return true;
    }
    // Get the optimal work size.
    worksize = -1;
    dgetri_(&n, V, &n, ipiv, &workQuery, &worksize, &info);
    CxmAssert(info == 0);
    worksize = workQuery;
    {
	double work[worksize];

	dgetri_(&n, V, &n, ipiv, work, &worksize, &info);
	if (info != 0) {
	    return true;
	}
    }
    for (int i = 0; i < n; i++) {
	for (int j = 0; j < n; j++) {
	    for (int k = 0; k < n; k++) {
		qEigVecCube[i*nSq + j*n + k] = QEigVecs[i + k*n] * V[k + j*n];
	    }
	}
    }

#ifdef CxmLikDebug
    fprintf(stderr, "VQ(V^-1):\n");
    for (unsigned i = 0; i < n; i++) {
	for (unsigned j = 0; j < n; j++) {
	    double p = 0.0;
	    for (int k = 0; k < n; k++) {
		p += qEigVecCube[i*nSq + j*n + k] * qEigVals[k];
	    }
	    fprintf(stderr, " %10.6f", p);
	}
	fprintf(stderr, "\n");
    }
#endif

    return false;
}

// Compute substitution probabilities, based on the eigen decomposition of a Q
// matrix, and branch length.  (Mutation rate is normalized in CxLikQ(), so
// that branch length has a standardized interpretation.)
//
//             Q*v
//   P(Q,v) = e
void
CxLikPt(int n, double *P, double *qEigVecCube, double *qEigVals, double v) {
    int nSq = n*n;
    double qEigValsExp[n];
    double p;

    for (int i = 0; i < n; i++) {
	qEigValsExp[i] = exp(qEigVals[i] * v);
    }

    for (int i = 0; i < n; i++) {
	for (int j = 0; j < n; j++) {
	    p = 0.0;
	    for (int k = 0; k < n; k++) {
		p += qEigVecCube[i*nSq + j*n + k] * qEigValsExp[k];
	    }
	    // Make sure there are no negative probabilities due to rounding
	    // error in eigen decomposition.
	    if (p < 0.0) {
		p = 0.0;
	    }
	    P[i*n + j] = p;
	}
    }

#ifdef CxmLikDebug
    fprintf(stderr, "----------------------------------------"
      "----------------------------------------\n");
    fprintf(stderr, "P(%f):\n", v);
    for (unsigned i = 0; i < n; i++) {
	for (unsigned j = 0; j < n; j++) {
	    fprintf(stderr, " %12.6e", P[i*n + j]);
	}
	fprintf(stderr, "\n");
    }
#endif
}

CxmpInline void
CxLikExecuteStripe(CxtLik *lik, unsigned stripe) {
    unsigned dim = lik->dim;
    unsigned dimSq = dim * dim;
    unsigned ncat = lik->ncat;
    unsigned dn = dim * ncat;
    unsigned cMin = lik->stripeWidth * stripe;
    unsigned cLim = cMin + lik->stripeWidth;
    double P[ncat][dimSq];

    // Iteratively process the execution plan.
    for (unsigned s = 0; s < lik->stepsLen; s++) {
	CxtLikStep *step = &lik->steps[s];
	double *parentMat = step->parentCL->cLMat;
	double *childMat = step->childCL->cLMat;
	unsigned glen = step->model->glen;

	// Compute P matrices, one for each gamma-distributed rate.
	for (unsigned g = 0; g < glen; g++) {
	    CxLikPt(dim, P[g], step->model->qEigVecCube, step->model->qEigVals,
	      step->edgeLen * step->model->gammas[g] * step->model->rmult *
	      step->model->wNorm);
	}

	if (step->variant == CxeLikStepComputeCL) {
	    for (unsigned c = cMin; c < cLim; c++) {
		for (unsigned g = 0; g < glen; g++) {
		    for (unsigned iP = 0; iP < dim; iP++) {
			double cL = 0.0;
			for (unsigned iC = 0; iC < dim; iC++) {
			    cL += P[g][iP*dim + iC] * childMat[c*dn + g*dim
			      + iC];
			}
			parentMat[c*dn + g*dim + iP] = cL;
		    }
		}
	    }
	} else {
	    CxmAssert(step->variant == CxeLikStepMergeCL);
	    for (unsigned c = cMin; c < cLim; c++) {
		for (unsigned g = 0; g < glen; g++) {
		    for (unsigned iP = 0; iP < dim; iP++) {
			double cL = 0.0;
			for (unsigned iC = 0; iC < dim; iC++) {
			    cL += P[g][iP*dim + iC] * childMat[c*dn + g*dim
			      + iC];
			}
			parentMat[c*dn + g*dim + iP] *= cL;
		    }
		}
	    }
	}
    }

    // For each model in the mixture, aggregrate conditional likelihoods and
    // weight them according to frequency priors, in order to compute
    // site-specific likelihoods.
    for (unsigned m = 0; m < lik->modelsLen; m++) {
	CxtLikModel *model = &lik->models[m];

	if (model->weightScaled == 0.0 || model->entire) {
	    // Model is not part of the mixture, or its siteL vector is
	    // already computed.
	    continue;
	}

	// Pre-compute state-specific weights, taking into account discretized
	// gamma-distributed rates.
	unsigned glen = model->glen;
	double pn[dim];
	for (unsigned i = 0; i < dim; i++) {
	    pn[i] = model->piDiagNorm[i] / (double)glen;
	}

	for (unsigned c = cMin; c < cLim; c++) {
	    double L = 0.0;
	    for (unsigned g = 0; g < glen; g++) {
		for (unsigned i = 0; i < dim; i++) {
		    L += pn[i] * model->cL->cLMat[c*dn + g*dim + i];
		}
	    }
	    model->siteL[c] = L;
	}
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
	CxMqPut(&CxpLikDoneMq, msg);
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
    if (CxMqNew(&CxpLikTodoMq, sizeof(CxtLikMsg *), CxNcpus * CxmLikMqMult)) {
	return;
    }
    if (CxMqNew(&CxpLikDoneMq, sizeof(CxtLikMsg *), CxNcpus * CxmLikMqMult)) {
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

CxmpInline void
CxpLikExecute(CxtLik *lik) {
    if (lik->stepsLen > 0) {
	if (CxNcpus > 1) {
	    pthread_once(&CxpLikOnce, CxpLikThreaded);
	}

	// Compute log-likelihoods for each model in the mixture, stripewise.
	if (CxpLikNThreads > 0 && lik->nstripes > 1) {
	    CxtLikMsg msgs[CxNcpus * CxmLikMqMult];
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
		CxmAssert(msg->lik == lik);
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
}

double
CxLikExecuteLnL(CxtLik *lik) {
    CxpLikExecute(lik);

    // Compute the weighted log-likelihood.
    double lnL = 0.0;
    for (unsigned c = 0; c < lik->nchars - lik->npad; c++) {
	double L = 0.0;
	for (unsigned m = 0; m < lik->modelsLen; m++) {
	    CxtLikModel *model = &lik->models[m];
	    if (model->weightScaled == 0.0) {
		// Model is not part of the mixture.
		continue;
	    }
	    L += model->siteL[c] * model->weightScaled;
	}
	lnL += log(L) * (double)lik->charFreqs[c];
    }
    CxmAssert(isnan(lnL) == 0);

    return lnL;
}

void
CxLikExecuteSiteLnLs(CxtLik *lik, double *lnLs) {
    CxpLikExecute(lik);

    for (unsigned c = 0; c < lik->nchars - lik->npad; c++) {
	double L = 0.0;
	for (unsigned m = 0; m < lik->modelsLen; m++) {
	    CxtLikModel *model = &lik->models[m];
	    if (model->weightScaled == 0.0) {
		// Model is not part of the mixture.
		continue;
	    }
	    L += model->siteL[c] * model->weightScaled;
	}
	double lnL = log(L);
	CxmAssert(isnan(lnL) == 0);
	lnLs[c] = lnL * (double)lik->charFreqs[c];
    }
}
