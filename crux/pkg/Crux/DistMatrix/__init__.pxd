# Forward declaration.
cdef class DistMatrix

from Crux.Tree cimport Tree
cimport Crux.Taxa as Taxa

from libc cimport *

cdef class DistMatrix:
    cdef readonly Taxa.Map taxaMap
    cdef float *dists
    cdef readonly size_t ntaxa
    cdef readonly bint additive # Unspecified unless rnj() has been used.

    cdef void _parse(self, input) except *
    cdef void _allocDists(self, size_t ntaxa) except *
    cdef void _dup(self, DistMatrix other, int sampleSize) except *
    cpdef float distanceGet(self, size_t x, size_t y)
    cpdef distanceSet(self, size_t x, size_t y, float distance)
    cdef Tree _nj(self, bint random)
    cdef Tree _rnj(self, bint random, bint additive)
    cdef void _rowsSwap(self, size_t a, size_t b)
    cdef void _matrixShuffle(self, list order) except *
    cpdef shuffle(self)
    cpdef Tree nj(self, bint joinRandom=*, bint destructive=*)
    cpdef Tree rnj(self, bint joinRandom=*, bint tryAdditive=*,
      bint destructive=*)
    cpdef render(self, str format=*, str distFormat=*, file outFile=*)

cdef inline size_t nxy2i(size_t n, size_t x, size_t y):
    assert x < n
    assert y < n
    assert x != y

    if x > y:
        x, y = y, x

    return n*x + y - (((x+3)*x) >> 1) - 1
