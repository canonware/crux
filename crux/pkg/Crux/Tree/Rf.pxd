from Crux.Tree cimport Tree, Ring

cdef class Vec:
    cdef unsigned nBits
    cdef unsigned char *bits

    cdef int cmp(self, Vec other)
    cdef reset(self)
    cdef bint get(self, unsigned bit)
    cdef void set(self, unsigned bit, bint val)
    cdef void invert(self)
    cdef void merge(self, Vec other)

cdef class Rf:
    cdef dict taxa
    cdef list edgeVecs
    cdef Vec leafVec

    cdef Vec _bipartitionsRecurse(self, Ring ring, bint calcVec)
    cdef void _bipartitions(self, Tree tree) except *
    cdef double dist(self, Rf other) except -1.0

cdef double rf(Tree a, Tree b) except -1.0
cdef list rfs(Tree a, list others)