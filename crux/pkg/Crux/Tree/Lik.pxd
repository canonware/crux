from Crux.Character cimport Character

cdef class QMat:
    cdef bint up2date
    cdef Character char_
    cdef unsigned dim
    cdef double *rMat
    cdef double *piMat
    cdef double *qMat
    cdef double *eigVecCube
    cdef double *eigVals

    cdef double getRate(self, unsigned i, unsigned j) except -1.0
    cdef void setRate(self, unsigned i, unsigned j, double rate) except *

    cdef double getFreq(self, unsigned i) except -1.0
    cdef void setFreq(self, unsigned i, double freq) except *

    cdef void _update(self) except *

cdef class PMat:
    cdef Character char_
    cdef unsigned dim
    cdef double *pMat
    cdef double *eigValsExp

    cdef void compute(self, QMat Q, double mu, double t) except *
