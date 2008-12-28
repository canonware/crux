from Crux.Tree cimport Tree, Node
cimport Crux.Taxa as Taxa

from libc cimport *

cdef class Nj:
    cdef float *dBase, *d # d is advanced as rows are removed.
    cdef size_t nBase, n
    cdef float *rBase, *r
    cdef float *rScaledBase, *rScaled
    cdef Tree tree
    cdef list nodes

#    cdef void _njDump(self) except *
    cdef void _rInit(self) except *
    cdef void _rScaledInit(self) except *
    cdef void _nodesInit(self, Taxa.Map taxaMap) except *
    cdef void _rScaledUpdate(self)
    cdef void _njRandomMinFind(self, size_t *rXMin, size_t *rYMin)
    cdef void _njDeterministicMinFind(self, size_t *rXMin, size_t *rYMin)
    cdef Node _njNodesJoin(self, size_t xMin, size_t yMin, float *rDistX,
      float *rDistY)
    cdef void _njRSubtract(self, size_t xMin, size_t yMin)
    cdef void _njCompact(self, size_t xMin, size_t yMin, Node node, float distX,
      float distY) except *
    cdef void _njDiscard(self)
    cdef void _njFinalJoin(self) except *

    cdef void prepare(self, float *d, size_t n, Taxa.Map taxaMap) except *
    cdef Tree nj(self, bint random)

cdef class Rnj(Nj):
    cdef size_t _rnjRowAllMinFind(self, size_t x, float *rDist)
    cdef bint _rnjRowAllMinOk(self, size_t x, float minDist)
    cdef size_t _rnjRowMinFind(self, size_t x)
    cdef bint _rnjPairClusterAdditive(self, size_t a, size_t b)
    cdef bint _rnjPairClusterOk(self, size_t a, size_t b)
    cdef bint _rnjRandomCluster(self, bint additive) except -1
    cdef bint _rnjDeterministicCluster(self, bint additive) except -1
    cdef rnj(self, bint random, bint additive)
