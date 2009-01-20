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
    // assuming DNA:
    //
    //     A C G T
    //   -----------
    //   | x x x x | 0
    //   | x x x x | 1
    //   | x x x x | 2
    //   | ....... | ...
    //   | x x x x | n-1
    //   -----------
    double *mat;
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

    // R, and Pi matrices, where Q = R*Pi (with diagonal set such that each row
    // sums to 0).
    double *rMat;
    double *piMat;

    // Q matrix eigenvector/eigenvalue decomposition.
    double *qEigVecCube;
    double *qEigVals;

    // Gamma-distributed mutation rates parameter, with discretization.  If
    // alpha is INFINITY, variable rates are disabled for this model.
    double alpha;
    double *gammas;

    // Site-specific conditional likelihoods for the root node.
    CxtLikCL cL;

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
    double *parentMat;
    double *childMat;
    double edgeLen;
} CxtLikStep;

typedef struct {
    // Number of character states (i.e. dimensionality of Q matrix).
    unsigned dim;

    // Number of discrete Gamma-distributed rate categories.
    unsigned ncat;

    // Number of characters.
    unsigned nchars;

    // Stripe width and number of stripes (stripe*nstripes == nchars).
    unsigned stripe;
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
    // plan.  At most, we require ((2*ntaxa - 3) * nmodels) steps (one per
    // branch per model).  stepsLen records the length of the current update
    // plan, and stepsMax records how much space is available.
    CxtLikStep *steps;
    unsigned stepsLen;
    unsigned stepsMax;
} CxtLik;

bool
CxLikQDecomp(int n, double *R, double *Pi, double *qEigVecCube,
  double *qEigVals);
void
CxLikPt(int n, double *P, double *qEigVecCube, double *qEigVals, double muT);

#endif // CxLik_h
