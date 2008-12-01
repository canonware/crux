# Forward declarations.
cdef class Tree
cdef class Node
cdef class Edge
cdef class Ring

from Crux.CTMatrix cimport CTMatrix
from Crux.Taxa cimport Taxon
cimport Crux.Taxa as Taxa

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
    # property ntaxa
    # property nnodes
    # property nedges
    # property base

    cpdef deroot(self)
    cpdef canonize(self, Taxa.Map taxaMap)
    cpdef int collapse(self) except -1
    cpdef tbr(self, Edge bisect, Edge reconnectA, Edge reconnectB)
    cpdef int tbrNNeigbhorsGet(self)
    cpdef tbrNeighborGet(self, int neighbor)
    cpdef nni(self, Edge edge, Edge reconnectA, Edge reconnectB)
    cpdef nniNNeigbhorsGet(self)
    cpdef nniNeighborGet(self, int neighbor)
    cpdef mpPrepare(self, CTMatrix cTMatrix, bint elimUninformative=*)
    cpdef mpFinish(self)
    cpdef mp(self)
    cpdef tbrBestNeighbhorsMp(self, int maxHold=*)
    cpdef tbrBetterNeighborsMp(self, int maxHold=*)
    cpdef tbrAllNeighborsMp(self, int maxHold=*)
    cpdef nHeldGet(self)
    cpdef heldGet(self, int i)
    cpdef str render(self, bint lengths=*, lengthFormat=*, Taxa.Map taxaMap=*)

cdef class Node:
    cdef object __weakref__
    cdef Tree _tree
    cdef Ring _ring
    cdef Taxon _taxon
    cdef int _degree

    # property tree
    # property taxon
    # property ring
    cpdef int _degreeGet(self, bint calculate=*) except -1
    # property degree
    cpdef rrender(self, Node prev, bint lengths, lengthFormat, Taxa.Map taxaMap,
      bint zeroLength=*, bint noLength=*)
    cpdef int separation(self, Node other)

cdef class Edge:
    cdef object __weakref__
    cdef Tree _tree
    cdef float _length
    cdef Ring _ring

    # property tree
    # property ring
    # property length
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
    cdef int _separation(self, Node other, int sep)

    # property tree
    # property node
    # property edge
    # property other
    # property next
    # property prev
