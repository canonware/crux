from libc cimport uint64_t

cdef extern from "CxLik.h":
    ctypedef struct CxtLikCL:
        uint64_t sn
        double *mat
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
        double *parentMat
        double *childMat
        double edgeLen
    ctypedef struct CxtLik:
        unsigned dim
        unsigned ncat
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

    cdef bint CxLikQDecomp(int n, double *R, double *Pi, double *qEigVecCube,
      double *qEigVals)
    cdef void CxLikPt(int n, double *P, double *qEigVecCube, double *qEigVals,
      double muRT)
    cdef bint CxLikExecute(CxtLik *lik, double *rLnL)
