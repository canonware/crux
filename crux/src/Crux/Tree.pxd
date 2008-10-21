# Forward declarations.
cdef class Tree
cdef class Node
cdef class Edge
cdef class Ring

from CTMatrix cimport CTMatrix

cdef class Tree:
    cdef dup(self)
    cdef taxonMapGet(self)
    cdef rf(self, Tree other)
    cdef ntaxaGet(self)
    cdef nedgesGet(self)
    cdef baseGet(self)
    cdef baseSet(self, Node base)
    cdef canonize(self)
    cdef collapse(self)
    cdef tbr(self, Edge bisect, Edge reconnectA, Edge reconnectB)
    cdef int tbrNNeigbhorsGet(self)
    cdef tbrNeighborGet(self, int neighbor)
    cdef nni(self, Edge edge, Edge reconnectA, Edge reconnectB)
    cdef nniNNeigbhorsGet(self)
    cdef nniNeighborGet(self, int neighbor)
    cdef mpPrepare(self, CTMatrix cTMatrix, bint elimUninformative=?)
    cdef mpFinish(self)
    cdef mp(self)
    cdef tbrBestNeighbhorsMp(self, int maxHold=?)
    cdef tbrBetterNeighborsMp(self, int maxHold=?)
    cdef tbrAllNeighborsMp(self, int maxHold=?)
    cdef nHeldGet(self)
    cdef heldGet(self, int i)
    cdef render(self, bint rooted=?, bint labels=?, bint lengths=?,
      lengthFormat=?, outFile=?)

cdef class Node:
    cdef Tree tree(self)
    cdef int taxonNumGet(self)
    cdef void taxonNumSet(self, int taxonNum)
    cdef Ring ring(self)
    cdef int degree(self)

cdef class Edge:
    cdef Tree _tree
    cdef float _length
    cdef Ring _ringA, _ringB

    cdef Tree tree(self)
    cdef rings(self)
    cdef float lengthGet(self)
    cdef void lengthSet(self, float length)
    cdef void attach(self, Node nodeA, Node nodeB)
    cdef detach(self)

cdef class Ring:
    cdef Node _node
    cdef Edge _edge
    cdef Ring _next, _prev

    cdef Tree tree(self)
    cdef Node node(self)
    cdef Edge edge(self)
    cdef Ring other(self)
    cdef Ring next(self)
    cdef Ring prev(self)
