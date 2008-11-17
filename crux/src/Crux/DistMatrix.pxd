from Tree cimport Tree
cimport Taxa

cdef class DistMatrix:
    cdef readonly Taxa.Map taxaMap
    cdef float *dists


    cdef void _parse(self, input) except *
    # property ntaxa
    # property taxaMap
    cdef void _dup(self, DistMatrix other, int sampleSize) except *
    cdef float distanceGet(self, int x, int y)
    cdef void distanceSet(self, int x, int y, float distance)
    cdef void _matrixShuffle(self, list order)
    cdef Tree _nj(self, bint random)
    cdef Tree _rnj(self, bint random, bint additive)
    cdef _render(self, format, distFormat, file file_=?)
