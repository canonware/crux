from libc cimport uint64_t

cdef extern from "CxLik.h":
    # Forward declarations.
    ctypedef struct CxtLikCL
    ctypedef struct CxtLikModel

    ctypedef struct CxtLikCL:
        double *cLMat
        double *lnScale
        bint valid
        CxtLikCL *parent
        unsigned nSibs
        double edgeLen
    ctypedef struct CxtLikComp:
        CxtLikModel *model
        double weightScaled
        double cweight
        double cmult
    ctypedef struct CxtLikModel:
        bint decomp
        double weight
        double qNorm
        double rmult
        unsigned *rclass
        double *rTri
        double *piDiag
        double *piDiagNorm
        double *qEigVecCube
        double *qEigVals
        double alpha
        bint catMedian
        bint invar
        unsigned comp0
        unsigned clen
    ctypedef enum CxeLikStep:
        CxeLikStepComputeL = 0
        CxeLikStepComputeI = 1
        CxeLikStepMergeL   = 2
        CxeLikStepMergeI   = 3
    ctypedef struct CxtLikStep:
        CxeLikStep variant
        unsigned ntrail
        CxtLikCL *parentCL
        CxtLikCL *childCL
        double edgeLen
    ctypedef struct CxtLik:
        unsigned polarity
        unsigned dim
        unsigned rlen
        unsigned nchars
        unsigned npad
        unsigned *charFreqs
        unsigned stripeWidth
        unsigned nstripes
        bint invalidate
        bint resize
        bint reweight
        double wNorm
        CxtLikModel **models
        unsigned modelsLen
        unsigned modelsMax
        CxtLikComp *comps
        unsigned compsLen
        unsigned compsMax
        CxtLikCL *rootCLC
        double *siteLnL
        CxtLikStep *steps
        unsigned stepsLen
        unsigned stepsMax

    cdef unsigned CxmLikMqMult

    cdef bint CxLikQDecomp(int n, double *RTri, double *PiDiag, \
      double *PiDiagNorm, double *qEigVecCube, double *qEigVals, double *qNorm)
    cdef void CxLikPt(int n, double *P, double *qEigVecCube, double *qEigVals, \
      double v)
    cdef void CxLikExecute(CxtLik *lik)
