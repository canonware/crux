from Tree cimport Tree
from TaxonMap cimport TaxonMap

cdef class DistMatrix:
    cdef void _parse(self, input, TaxonMap taxonMap, bint symmetric)
    cdef int ntaxaGet(self)
    cdef bint isSymmetric(self)
    cdef DistMatrix _dup(self, input, TaxonMap taxonMap)
    cdef DistMatrix _sample(self, input, TaxonMap taxonMap, list rows)
    cdef TaxonMap taxonMapGet(self)
    cdef float distanceGet(self, int x, int y)
    cdef void distanceSet(self, int x, int y, float distance)
    cdef void _matrixShuffle(self, list order)
    cdef Tree _nj(self, bint random)
    cdef Tree _rnj(self, bint random, bint additive)
    cdef _render(self, format, distFormat, file file_=?)
