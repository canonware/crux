# Forward declarations.
cdef class Tree
cdef class Node
cdef class Edge
cdef class Ring

from CTMatrix cimport CTMatrix
from Taxa cimport Taxon
cimport Taxa

cdef class Tree:
    cdef Node _base
    cdef _taxa, _nodes, _edges
    cdef list _renderList
    cdef int _sn # Incremented every time the tree is modified.
    cdef int _cacheSn, _cachedNtaxa, _cachedNnodes, _cachedNedges
    cdef public bint rooted

    cdef void _randomNew(self, int ntaxa, Taxa.Map taxaMap) except *
    cdef void _newickNew(self, str input, Taxa.Map taxaMap) except *

    cdef Node _dup(self, Tree newTree, Node node, Ring prevRing)
    cpdef Tree dup(self)
    # property taxa
    # property nodes
    # property edges
    cpdef rf(self, Tree other)
    cdef void _recacheRecurse(self, Ring ring)
    cdef void _recache(self)
    cpdef int ntaxaGet(self)
    cpdef int nnodesGet(self)
    cpdef int nedgesGet(self)
    cpdef Node baseGet(self)
    cpdef baseSet(self, Node base)

    cpdef deroot(self)
    cpdef canonize(self, Taxa.Map taxaMap)
    cpdef int collapse(self) except -1
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
    cpdef str render(self, bint lengths=?, lengthFormat=?, Taxa.Map taxaMap=?)

cdef class Node:
    cdef object __weakref__
    cdef Tree _tree
    cdef Ring _ring
    cdef Taxon _taxon
    cdef int _degree

    cpdef Tree tree(self)
    # property taxon
    cpdef Ring ring(self)
    cpdef int degree(self, bint calculate=?) except -1
    cpdef rrender(self, Node prev, bint lengths, lengthFormat, Taxa.Map taxaMap,
      bint zeroLength=?, bint noLength=?)
    cpdef int separation(self, Node other)

cdef class Edge:
    cdef object __weakref__
    cdef Tree _tree
    cdef float _length
    cdef Ring _ring

    cpdef Tree tree(self)
    cpdef rings(self)
    cpdef float lengthGet(self)
    cpdef lengthSet(self, float length)
    cpdef attach(self, Node nodeA, Node nodeB)
    cpdef detach(self)

cdef class Ring:
    cdef Node _node
    cdef Edge _edge
    cdef Ring _other
    cdef Ring _next, _prev

    cpdef siblings(self) # Iterator.

    cdef Node _minTaxon(self, Taxa.Map taxaMap)
    cdef Node _canonize(self, Taxa.Map taxaMap)
    cdef void _collapsable(self, list collapsable) except *
    cdef void _collapse(self)

    cpdef Tree tree(self)
    cpdef Node node(self)
    cpdef Edge edge(self)
    cpdef Ring other(self)
    cpdef Ring next(self)
    cpdef Ring prev(self)
    cdef int _separation(self, Node other, int sep)
