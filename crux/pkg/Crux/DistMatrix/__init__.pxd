# Forward declaration.
cdef class DistMatrix

from Crux.Tree cimport Tree
cimport Crux.Taxa as Taxa

from CxDistMatrix cimport *

cdef class DistMatrix:
    cdef readonly Taxa.Map taxaMap
    cdef CxtDMDist *dists
    cdef readonly CxtDMSize ntaxa
    cdef readonly bint additive # Unspecified unless rnj() has been used.

    cdef void _parse(self, input) except *
    cdef void _allocDists(self, CxtDMSize ntaxa) except *
    cdef void _dup(self, DistMatrix other, int sampleSize) except *
    cpdef CxtDMDist distanceGet(self, CxtDMSize x, CxtDMSize y)
    cpdef distanceSet(self, CxtDMSize x, CxtDMSize y, CxtDMDist distance)
    cdef Tree _nj(self, bint random)
    cdef Tree _rnj(self, bint random, bint additive)
    cdef void _rowsSwap(self, CxtDMSize a, CxtDMSize b)
    cdef void _matrixShuffle(self, list order)
    cpdef shuffle(self)
    cpdef Tree nj(self, bint joinRandom=*, bint destructive=*)
    cpdef Tree rnj(self, bint joinRandom=*, bint tryAdditive=*,
      bint destructive=*)
    cpdef render(self, str format=*, str distFormat=*, file outFile=*)
