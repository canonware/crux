#ifndef CxLik_h
#define CxLik_h

#include "../Cx.h"

typedef struct CxsLikCL CxtLikCL;
typedef struct CxsLikModel CxtLikModel;

struct CxsLikCL {
    // Conditional likelihood matrix, stored in row-major form, such that each
    // row contains conditional likelihoods for a single site.  For example,
    // assuming DNA and 4 model components (whether due to mixture models
    // and/or +G rate categories), cLMat for an internal node looks like:
    //
    //     comp 0    comp 1    comp 2    comp 3
    //     A C G T   A C G T   A C G T   A C G T
    //   -----------------------------------------
    //   | x x x x | x x x x | x x x x | x x x x | 0
    //   | x x x x | x x x x | x x x x | x x x x | 1
    //   | ....... | ....... | ....... | ....... | ...
    //   | x x x x | x x x x | x x x x | x x x x | n-1
    //   -----------------------------------------
    //
    // cLMat may be deallocated (and the pointer set to NULL) if execution
    // planning finds it to be obsolete.
    //
    // Leaf nodes only store one copy of the character data, regardless of the
    // number of model components:
    //
    //     A C G T
    //   -----------
    //   | x x x x | 0
    //   | x x x x | 1
    //   | ....... | ...
    //   | x x x x | n-1
    //   -----------
    double *cLMat;

    // Vector of character-specific log-scale factors.  The conditional
    // likelihoods for each character are rescaled such that the largest
    // conditional likelihood is 1.0.  This avoids floating point underflow
    // issues, but the rescaling has to be tracked so that it can be accounted
    // for in final lnL computation.
    double *lnScale;

    // True if the contents of cLMat and lnScale are consistent with the
    // current tree topology.  This field is cleared during recursive execution
    // planning if any child determines the topology is incompatible with its
    // parent's cache.
    bool valid;

    // Parent which most recently used this cache when computing its conditional
    // likelihood.  All siblings must still refer to the parent, the number of
    // siblings must still be nSibs, and the branches separating them from the
    // parent must remain the same length, in order for the parent's cache to be
    // potentially valid.  (The model of evolution must also remain unchanged.)
    CxtLikCL *parent;
    unsigned nSibs;
    double edgeLen;
};

// Model component.  Each model in the mixture consists of one or more
// model components.
typedef struct {
    // Model this component is associated with.
    CxtLikModel *model;

    // Scaled weight, such that weightScaled for all model components sums to 1.
    double weightScaled;

    // Component relative weight, with relation to other components associated
    // with the same model.
    double cweight;

    // Component rate multiplier.  This is 1 unless:
    //
    // * used for a +G rate category, in which case the mean cmult across the
    //   model's +G rate category components is 1, or
    // * used for a +I component, in which case it is 0.
    double cmult;
} CxtLikComp;

// One or more mixture models are used to compute a tree's lnL.  Each model is
// given a weight (scaled model component weights sum to 1) and a Q
// normalization scaler to adjust the mean rate to 1 across the entire mixture.
// Additionally, each model may have a discretized Gamma distribution of
// mutation rates (average rate of 1).  The precomputed eigenvector/eigenvalue
// decompostion of the Q matrix is used to compute the probability of
// substitutions along a branch of a specific length.
struct CxsLikModel {
    // If any model parameters have changed since the last time eigen
    // decomposition was performed, decomp is true.
    bool decomp;

    // Mixture model relative weight.
    double weight;

    // Where P(t) is the mutation probability matrix for time interval t,
    // (t == v*wNorm) (wNorm is in CxtLik), where v is branch length.  qNorm is
    // a normalization factor that scales qNorm*Q such that the mean
    // substitution rate is 1.  If only a single model is in use, (wNorm ==
    // qNorm); otherwise wNorm may deviate from qNorm in order to adjust the
    // weighted mean substitution rate to 1 across the entire mixture.
    double qNorm;

    // Upper triangle of the R matrix and the diagonal of the Pi matrix, where
    // Q = (rmult*R)*Pi (with diagonal set such that each row sums to 0).
    double rmult;
    unsigned *rclass;
    double *rTri;
    double *piDiag;
    double *piDiagNorm;

    // Q matrix eigenvector/eigenvalue decomposition.
    double *qEigVecCube;
    double *qEigVals;

    // Gamma-distributed mutation rates parameter.  If alpha is INFINITY,
    // variable rates effectively disabled for this model, and only the first
    // element of comps has a non-zero scaled weight.  Otherwise, each element
    // of comps has equal weight, and cmult is set appropriately in each
    // element to represent +G rates.
    double alpha;
    // Use category means if catMedian is false, category medians if true.
    bool catMedian;

    // If true, the model incorporates invariable sites via an extra model
    // component that is at index comp0+clen-1.  The proportion of invariable
    // sites is implicit in the cweight ratio between the invar component and
    // all the other components for this model.
    bool invar;

    // comp0 is the index of the first component within CxtLik's comps vector
    // that corresponds to this model; components for each model are contiguous
    // within the comps vector.  clen is the number of components used by this
    // model.
    unsigned comp0;
    unsigned clen;
};

// Updating the cache and computing the lnL of the tree can be broken down into
// two basic types of computations.  The conditional likelihood of a node is
// updated by computing the conditional likelihood of one child
// (CxeLikStepCompute*), then computing and merging the conditional
// likelihoods of other children, one at a time (CxeLikStepMerge*).
//
// Since CL's associated with leaves only store a single copy of the character
// data, "L" variants of the algorithm are necessary to handle those CL's,
// whereas "I" variants handle CL's associated with internal nodes.
typedef enum {
    CxeLikStepComputeL = 0,
    CxeLikStepComputeI = 1,
    CxeLikStepMergeL   = 2,
    CxeLikStepMergeI   = 3
} CxeLikStep;

// Each update step needs the same basic pieces of information.  The data that
// vary among steps are stored using this data structure.  parentCL and childCL
// are the conditional likelihood data associated with the node being updated,
// and one of its children.  edgeLen is the length of the branch separating the
// two nodes.
//
//        node
//         |
//         |edge
//         |
//       child
typedef struct {
    CxeLikStep variant;
    unsigned ntrail; // Number of trailing merge steps.
    CxtLikCL *parentCL;
    CxtLikCL *childCL;
    double edgeLen;
} CxtLikStep;

typedef struct {
    // Polarity of data structures (0 or 1).  It is possible for a Lik to have
    // a clone (created via Lik.clone()) that uses the opposite polarity.
    unsigned polarity;

    // Number of character states (i.e. dimensionality of Q matrix).
    unsigned dim;

    // Maximum number of independent relative mutation rate classes (i.e.
    // length of rclass and rTri).
    unsigned rlen;

    // Number of characters.
    unsigned nchars;

    // Number of pad characters (used to make nchars a multiple of stripeWidth).
    unsigned npad;

    // Character frequencies for the compact alignment.
    unsigned *charFreqs;

    // Stripe width and number of stripes (stripe*nstripes == nchars).
    unsigned stripeWidth;
    unsigned nstripes;

    // True if the cached CL's have been rendered invalid by some model change.
    // This flag is cleared during execution planning, which is when cache
    // invalidation happens.
    bool invalidate;

    // True if the mixture has grown/shrunk, which means that CL's need to
    // re-size their cLMat's.
    bool resize;

    // True if any relative weights have changed for the mixture models since
    // the last time wNorm was computed.
    bool reweight;

    // See CxtLikModel's qNorm for explanation.
    double wNorm;

    // Mixture models vector.  modelsLen indicates how many models are
    // currently in the mixture, and modelsMax indicates the total number of
    // model slots available without reallocating.
    CxtLikModel **models;
    unsigned modelsLen;
    unsigned modelsMax;

    // Array of model components correpsonding to the mixture models vector.
    // compsLen indicates how many components are currently in use, and
    // compsMax indicates the total number of components available in the
    // vector without reallocating.  Since models may use varying numbers of
    // model components, modelsLen is not necessarily the same as compsLen.
    CxtLikComp *comps;
    unsigned compsLen;
    unsigned compsMax;

    // Site-specific conditional likelihoods for the root node.
    CxtLikCL *rootCLC;

    // Array in which to place site likelihoods.
    double *siteLnL;

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

// Use message queues that have CxNcpus * CxmLikMqMult slots to communicate
// with worker threads.
#define CxmLikMqMult 8

bool
CxLikQDecomp(int n, double *RTri, double *PiDiag, double *PiDiagNorm,
  double *qEigVecCube, double *qEigVals, double *qNorm);
void
CxLikPt(int n, double *P, double *qEigVecCube, double *qEigVals, double v);
void
CxLikExecute(CxtLik *lik);

#endif // CxLik_h
