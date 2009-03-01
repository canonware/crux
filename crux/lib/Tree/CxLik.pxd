from libc cimport uint64_t

cdef extern from "CxLik.h":
    ctypedef struct CxtLikCL:
        uint64_t sn
        double *cLMat
    ctypedef struct CxtLikModel:
        uint64_t sn
        bint reassign
        double weight
        double weightScaled
        double qNorm
        double wNorm
        double rmult
        unsigned *rclass
        double *rTri
        double *piDiag
        double *piDiagNorm
        double *qEigVecCube
        double *qEigVals
        double alpha
        double *gammas
        unsigned glen
        CxtLikCL *cL
        bint entire
        double *siteL
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
        unsigned rlen
        unsigned ncat
        bint catMedian
        unsigned nchars
        unsigned npad
        unsigned *charFreqs
        unsigned stripeWidth
        unsigned nstripes
        bint renorm
        CxtLikModel *models
        unsigned modelsLen
        unsigned modelsMax
        CxtLikStep *steps
        unsigned stepsLen
        unsigned stepsMax

    cdef unsigned CxmLikMqMult

    cdef bint CxLikQDecomp(int n, double *RTri, double *PiDiag, \
      double *PiDiagNorm, double *qEigVecCube, double *qEigVals, double *qNorm)
    cdef void CxLikPt(int n, double *P, double *qEigVecCube, double *qEigVals, \
      double v)
    cdef double CxLikExecuteLnL(CxtLik *lik)
    cdef void CxLikExecuteSiteLnLs(CxtLik *lik, double *lnLs)
