# Forward declarations.
cdef class Tree
cdef class Node
cdef class Edge
cdef class Ring

from CTMatrix cimport CTMatrix
from TaxonMap cimport TaxonMap

cdef class Tree:
    cdef TaxonMap _taxonMap
    cdef _base

    cpdef dup(self)
    cpdef taxonMapGet(self)
    cpdef rf(self, Tree other)
    cpdef ntaxaGet(self)
    cpdef nedgesGet(self)
    cpdef baseGet(self)
    cpdef baseSet(self, Node base)
    cpdef canonize(self)
    cpdef collapse(self)
    cpdef tbr(self, Edge bisect, Edge reconnectA, Edge reconnectB)
    cpdef int tbrNNeigbhorsGet(self)
    cpdef tbrNeighborGet(self, int neighbor)
    cpdef nni(self, Edge edge, Edge reconnectA, Edge reconnectB)
    cpdef nniNNeigbhorsGet(self)
    cpdef nniNeighborGet(self, int neighbor)
    cpdef mpPrepare(self, CTMatrix cTMatrix, bint elimUninformative=?)
    cpdef mpFinish(self)
    cpdef mp(self)
    cpdef tbrBestNeighbhorsMp(self, int maxHold=?)
    cpdef tbrBetterNeighborsMp(self, int maxHold=?)
    cpdef tbrAllNeighborsMp(self, int maxHold=?)
    cpdef nHeldGet(self)
    cpdef heldGet(self, int i)
    cpdef render(self, bint rooted=?, bint labels=?, bint lengths=?,
      lengthFormat=?, outFile=?)

cdef class Node:
    cdef Tree _tree
    cdef Ring _ring
    cdef int _taxonNum

    cpdef Tree tree(self)
    cpdef int taxonNumGet(self)
    cpdef taxonNumSet(self, int taxonNum)
    cpdef Ring ring(self)
    cpdef int degree(self)

cdef class Edge:
    cdef Tree _tree
    cdef float _length
    cdef Ring _ringA, _ringB

    cpdef Tree tree(self)
    cpdef rings(self)
    cpdef float lengthGet(self)
    cpdef lengthSet(self, float length)
    cpdef attach(self, Node nodeA, Node nodeB)
    cpdef detach(self)

cdef class Ring:
    cdef Node _node
    cdef Edge _edge
    cdef Ring _next, _prev

    cpdef Tree tree(self)
    cpdef Node node(self)
    cpdef Edge edge(self)
    cpdef Ring other(self)
    cpdef Ring next(self)
    cpdef Ring prev(self)
