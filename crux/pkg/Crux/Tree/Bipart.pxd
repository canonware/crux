# Forward declarations.
cdef class Vec
cdef class Bipart

from Crux.Tree cimport Tree, Edge, Ring

cdef class Vec:
    cdef readonly Edge edge
    cdef readonly unsigned nBits
    cdef unsigned nBytes
    cdef unsigned char *bits

    cpdef int cmp(self, Vec other) except *
    cpdef reset(self)
    cpdef bint get(self, unsigned bit) except *
    cdef _set(self, unsigned bit, bint val)
    cpdef set(self, unsigned bit, bint val)
    cpdef invert(self)
    cpdef merge(self, Vec other)

cdef class Bipart:
    cdef dict taxaX # taxon-->index translation
    cdef readonly bint leaves # If True, edgeVecs includes leaf edges.
    cdef readonly list edgeVecs
    cdef Vec leafVec # Temp vector.

    cpdef int cmp(Bipart self, Bipart other)
    cdef Vec _bipartitionsRecurse(self, Ring ring, bint calcVec)
    cdef void _bipartitions(self, Tree tree) except *
    cpdef double rfDist(self, Bipart other) except -1.0
