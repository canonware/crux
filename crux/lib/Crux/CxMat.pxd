cdef extern from "CxMat.h":
    cdef inline double CxMatDdet(unsigned n, double *A)
    cdef inline double CxMatLogDet(unsigned n, double *A)
