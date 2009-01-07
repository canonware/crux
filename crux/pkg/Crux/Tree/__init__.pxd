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
    cdef file _renderFile
    cdef int _sn # Incremented every time the tree is modified.
    cdef int _cacheSn
    cdef list _cachedTaxa, _cachedNodes, _cachedEdges
    cdef public bint rooted

    cdef void _randomNew(self, int ntaxa, Taxa.Map taxaMap) except *
    cdef void _newickNew(self, str input, Taxa.Map taxaMap) except *

    cdef Node _dup(self, Tree newTree, Node node, Ring prevRing)
    cpdef Tree dup(self)
    cpdef double rf(self, Tree other)
    cpdef list rfs(self, list others)
    cdef void _recacheRecurse(self, Ring ring)
    cdef void _recache(self)
    cdef int getNtaxa(self)
    # property ntaxa
    cdef int getNnodes(self) except *
    # property nnodes
    cdef int getNedges(self) except *
    # property nedges
    cdef list getTaxa(self)
    # property taxa
    cdef list getNodes(self)
    # property nodes
    cdef list getEdges(self)
    # property edges
    cdef Node getBase(self)
    cdef void setBase(self, Node base)
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
    cdef void _renderAppend(self, str s) except *
    cpdef str render(self, bint lengths=*, lengthFormat=*, Taxa.Map taxaMap=*,
      file outFile=*)

cdef class Node:
    cdef object __weakref__
    cdef readonly Tree tree
    cdef readonly Ring ring
    cdef Taxon _taxon
    cdef int _degree

    cdef Taxon getTaxon(self)
    cdef void setTaxon(self, Taxon taxon) except *
    # property taxon
    cpdef int _degreeGet(self, bint calculate=*) except -1
    cdef int getDegree(self)
    # property degree
    cpdef rrender(self, Node prev, bint lengths, lengthFormat, Taxa.Map taxaMap,
      bint zeroLength=*, bint noLength=*)
    cpdef int separation(self, Node other)

cdef class Edge:
    cdef object __weakref__
    cdef readonly Tree tree
    cdef double _length
    cdef readonly Ring ring

    cdef double getLength(self)
    cdef void setLength(self, double length)
    # property length
    cpdef attach(self, Node nodeA, Node nodeB)
    cpdef detach(self)

cdef class Ring:
    cdef readonly Node node
    cdef readonly Edge edge
    cdef readonly Ring other
    cdef readonly Ring next
    cdef readonly Ring prev

    cpdef siblings(self) # Iterator.

    cdef Node _minTaxon(self, Taxa.Map taxaMap)
    cdef Node _someLeaf(self)
    cdef Node _canonize(self, Taxa.Map taxaMap)
    cdef void _collapsable(self, list collapsable, list clampable) except *
    cdef void _collapse(self)
    cdef int _separation(self, Node other, int sep)
