from Tree cimport Tree
cimport Taxa

from CxDistMatrix cimport *

cdef class DistMatrix:
    cdef readonly Taxa.Map taxaMap
    cdef CxtDMDist *dists
    cdef readonly CxtDMSize ntaxa

    cdef void _parse(self, input) except *
    cdef void _allocDists(self, CxtDMSize ntaxa) except *
    cdef void _dup(self, DistMatrix other, int sampleSize) except *
    cdef CxtDMDist distanceGet(self, CxtDMSize x, CxtDMSize y)
    cdef void distanceSet(self, CxtDMSize x, CxtDMSize y, CxtDMDist distance)
    cdef void _matrixShuffle(self, list order)
    cdef Tree _nj(self, bint random)
    cdef Tree _rnj(self, bint random, bint additive)
    cpdef render(self, str format=?, str distFormat=?, file outFile=?)
