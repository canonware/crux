cdef extern from "CxMat.h":
    cdef inline double CxMatDdet(unsigned n, double *A)
    cdef inline double CxMatLogDet(unsigned n, double *A)
    cdef void CxMatQDecomp(unsigned n, double *Q, double *eigVecCube, \
      double *eigVals)
    cdef void CxMatPt(int n, double *P, double *eigVecCube, double *eigVals, \
      double muT)
