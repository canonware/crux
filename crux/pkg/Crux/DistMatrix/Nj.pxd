from Crux.Tree cimport Tree, Node
cimport Crux.Taxa as Taxa

from CxDistMatrix cimport *
from CxDistMatrixNj cimport *

cdef class Nj:
    cdef CxtDMDist *dBase, *d # d is advanced as rows are removed.
    cdef CxtDMSize nBase, n
    cdef CxtDMDist *rBase, *r
    cdef CxtDMDist *rScaledBase, *rScaled
    cdef Tree tree
    cdef list nodes

    cdef void _rInit(self) except *
    cdef void _rScaledInit(self) except *
    cdef void _nodesInit(self, Taxa.Map taxaMap) except *
    cdef void _rScaledUpdate(self)
    cdef void _njRandomMinFind(self, CxtDMSize *rXMin, CxtDMSize *rYMin)
    cdef void _njDeterministicMinFind(self, CxtDMSize *rXMin, CxtDMSize *rYMin)
    cdef Node _njNodesJoin(self, CxtDMSize xMin, CxtDMSize yMin,
      CxtDMDist *rDistX, CxtDMDist *rDistY)
    cdef void _njRSubtract(self, CxtDMSize xMin, CxtDMSize yMin)
    cdef void _njCompact(self, CxtDMSize xMin, CxtDMSize yMin, Node node,
      CxtDMDist distX, CxtDMDist distY) except *
    cdef void _njDiscard(self)

    cdef void prepare(self, CxtDMDist *d, CxtDMSize n, Taxa.Map taxaMap) \
      except *
    cdef Tree nj(self, bint random)

cdef class Rnj(Nj):
    cdef Tree rnj(self, bint random, bint additive)
