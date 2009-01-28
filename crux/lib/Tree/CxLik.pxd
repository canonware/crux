from libc cimport uint64_t

cdef extern from "CxLik.h":
    ctypedef struct CxtLikCL:
        uint64_t sn
        double *cLMat
        double *lnScale
    ctypedef struct CxtLikModel:
        uint64_t sn
        bint reassign
        double weight
        double *rMat
        double *piDiag
        double *qEigVecCube
        double *qEigVals
        double alpha
        double *gammas
        unsigned glen
        CxtLikCL *cL
        bint entire
        double *stripeLnL
        double lnL
    ctypedef enum CxeLikStep:
        CxeLikStepComputeCL = 0
        CxeLikStepMergeCL   = 1
    ctypedef struct CxtLikStep:
        CxeLikStep variant
        CxtLikModel *model
        CxtLikCL *parentCL
        CxtLikCL *childCL
        double edgeLen
    ctypedef struct CxtLik:
        unsigned dim
        unsigned ncat
        bint catMedian
        unsigned nchars
        unsigned *charFreqs
        unsigned stripeWidth
        unsigned nstripes
        CxtLikModel *models
        unsigned modelsLen
        unsigned modelsMax
        CxtLikStep *steps
        unsigned stepsLen
        unsigned stepsMax

    cdef unsigned CxmLikMqMult

    cdef bint CxLikQDecomp(int n, double *R, double *Pi, double *qEigVecCube,
      double *qEigVals)
    cdef void CxLikPt(int n, double *P, double *qEigVecCube, double *qEigVals,
      double v)
    cdef double CxLikExecute(CxtLik *lik)
