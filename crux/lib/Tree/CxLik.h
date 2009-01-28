#ifndef CxLik_h
#define CxLik_h

#include "../Cx.h"

// Each ring object in a tree has a vector of CxtLikCL elements
// associated with it, one element for each mixture model.
typedef struct {
    // Model serial number associated with the cached conditional likelihood
    // matrix.  In order for the cache to be usable, the stored serial number
    // must match that of the model's current serial number.
    uint64_t sn;

    // Conditional likelihood matrix, stored in row-major form, such that each
    // row contains conditional likelihoods for a single site.  For example,
    // assuming DNA and no gamma-distributed rate categories:
    //
    //     A C G T
    //   -----------
    //   | x x x x | 0
    //   | x x x x | 1
    //   | x x x x | 2
    //   | ....... | ...
    //   | x x x x | n-1
    //   -----------
    //
    // If gamma-distributed rate categories are enabled, the matrix contains
    // multiple rows for each site (one per rate category).  For example, with
    // 4 rate categories:
    //
    //     A C G T
    //   -----------
    //   | x x x x | 0 (cat 0
    //   | x x x x | 0      1
    //   | x x x x | 0      2
    //   | x x x x | 0      3)
    //   | x x x x | 1 (cat 0
    //   | x x x x | 1      1
    //   | x x x x | 1      2
    //   | x x x x | 1      3)
    //   | ....... | ...
    //   | x x x x | n-1 (cat 0
    //   | x x x x | n-1      1
    //   | x x x x | n-1      2
    //   | x x x x | n-1      3)
    //   -----------
    double *cLMat;

    // Vector of character-specific log-scale factors.  The conditional
    // likelihoods for each character are rescaled such that the largest
    // conditional likelihood is 1.0.  This avoids floating point underflow
    // issues, but the rescaling has to be tracked so that it can be accounted
    // for in final lnL computation.
    double *lnScale;
} CxtLikCL;

// One or more mixture models are used to compute a tree's lnL.  Each model is
// given a weight (weights sum to 1).  Additionally, each model may have a
// discretized Gamma distribution of mutation rates (average rate of 1).  The
// precomputed eigenvector/eigenvalue decompostion of the Q matrix is used to
// compute the probability of substitutions along a branch of a specific length.
typedef struct {
    // Model serial number.  sn is unique to each proposed model, and is used
    // as the key to determine whether a cached conditional likelihood matrix
    // (as stored in CxtLikCL) is potentially usable.
    uint64_t sn;

    // If any model parameters have changed since the last time the serial
    // number was assigned, reassign is true.
    bool reassign;

    // Mixture model weight.  This is 1.0 if only a single model is in use.
    double weight;

    // R matrix and the diagonal of the Pi matrix, where Q = R*Pi (with
    // diagonal set such that each row sums to 0).
    double *rMat;
    double *piDiag;

    // Q matrix eigenvector/eigenvalue decomposition.
    double *qEigVecCube;
    double *qEigVals;

    // Gamma-distributed mutation rates parameter, with discretization.  If
    // alpha is INFINITY, variable rates are disabled for this model, glen is
    // 1, and gammas[0] is 1.0.  Otherwise, (glen == ncat) and gammas[...] is
    // set appropriately.
    double alpha;
    double *gammas;
    unsigned glen;

    // Site-specific conditional likelihoods for the root node.
    CxtLikCL *cL;

    // False if execution planning finds any places where conditional
    // likelihoods have to be recomputed.  If true, then the contents of
    // stripeLnL and lnL are valid, which means that absolutely no computation
    // is necessary to determine the lnL under this model.
    bool entire;

    // Array in which to place the sum of each stripe's log-likelihoods.
    double *stripeLnL;

    // Tree's unweighted log-likelihood.
    double lnL;
} CxtLikModel;

// Updating the cache and computing the lnL of the tree can be broken down into
// two types of computations.  The conditional likelihood of a node is updated
// by computing the conditional likelihood of one child (StepComputeCL), then
// computing and merging the conditional likelihoods of other children, one at
// a time (StepMergeCL).
typedef enum {
    CxeLikStepComputeCL = 0,
    CxeLikStepMergeCL   = 1
} CxeLikStep;

// Each update step needs the same basic pieces of information.  The data that
// vary among steps are stored using this data structure.  parentMat and
// childMat are the conditional likelihood matrices associated with the node
// being updated, and one of its children.  edgeLen is the length of the branch
// separating the two nodes.
//
//        node
//         |
//         |edge
//         |
//       child
typedef struct {
    CxeLikStep variant;
    CxtLikModel *model;
    CxtLikCL *parentCL;
    CxtLikCL *childCL;
    double edgeLen;
} CxtLikStep;

typedef struct {
    // Number of character states (i.e. dimensionality of Q matrix).
    unsigned dim;

    // Number of discrete Gamma-distributed rate categories.
    unsigned ncat;
    // Use category means if catMedian is false, category medians if true.
    bool catMedian;

    // Number of characters.
    unsigned nchars;

    // Character frequencies for the compact alignment.
    unsigned *charFreqs;

    // Stripe width and number of stripes (stripe*nstripes == nchars).
    unsigned stripeWidth;
    unsigned nstripes;

    // Mixture models vector.  modelsLen indicates how many models are
    // currently in the mixture, and modelsMax indicates the total number of
    // model slots available without reallocating.
    CxtLikModel *models;
    unsigned modelsLen;
    unsigned modelsMax;

    // Updating a tree's conditional likelihoods in order to compute the tree's
    // lnL is performed in two stages.  The first stage is a post-order tree
    // traversal that computes what steps must be taken in order to update the
    // tree's caches and complete the computation.  That information is stored
    // in this data structure, so that it is possible to (optionally) feed the
    // plan to a set of worker threads that perform the computations one stripe
    // at a time.
    //
    // The steps array size is incrementally increased as necessary so that
    // there is always enough space in the steps array to store the full update
    // plan.  At most, we require ((2*ntaxa - 2) * nmodels) steps (one per
    // branch per model, keeping in mind that we may need one extra per model
    // for the synthetic root).  stepsLen records the length of the current
    // update plan, and stepsMax records how much space is available.
    CxtLikStep *steps;
    unsigned stepsLen;
    unsigned stepsMax;
} CxtLik;

// Use message queues to communicate with worker threads that have CxNcpus *
// CxmLikMqMult slots.
#define CxmLikMqMult 8

bool
CxLikQDecomp(int n, double *R, double *Pi, double *qEigVecCube,
  double *qEigVals);
void
CxLikPt(int n, double *P, double *qEigVecCube, double *qEigVals, double v);
double
CxLikExecute(CxtLik *lik);

#endif // CxLik_h
