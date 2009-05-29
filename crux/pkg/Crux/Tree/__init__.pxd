# Forward declarations.
cdef class Tree
cdef class Node
cdef class Edge
cdef class Ring

from libc cimport *
from Crux.CTMatrix cimport CTMatrix
from Crux.Taxa cimport Taxon
cimport Crux.Taxa as Taxa
from Crux.Tree.Bipart cimport Bipart

cdef class Tree:
    cdef Node _base
    cdef list _renderList
    cdef file _renderFile
    # Incremented every time the tree is modified.
    cdef readonly uint64_t sn
    cdef int64_t _cacheSn
    cdef list _cachedTaxa, _cachedNodes, _cachedEdges
    cdef Bipart _cachedBipart
    cdef public bint rooted
    cdef public object aux

    cdef void _randomNew(self, int ntaxa, Taxa.Map taxaMap) except *
    cdef void _newickNew(self, str input, Taxa.Map taxaMap) except *

    cdef Node _dup(self, Tree newTree, Node node, Ring prevRing)
    cpdef Tree dup(self)
    cpdef double rf(self, Tree other) except -1.0
    cpdef list rfs(self, list others)
    cdef void _resetCache(self) except *
    cpdef clearCache(self)
    cdef void _recacheRecurse(self, Ring ring) except *
    cdef void _recache(self) except *
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
    cdef Bipart getBipart(self)
    # property bipart
    cdef Node getBase(self)
    cdef void setBase(self, Node base) except *
    # property base

    cdef void _clearAuxRecurse(self, Ring ring)
    cpdef clearAux(self)
    cpdef deroot(self)
    cpdef canonize(self, Taxa.Map taxaMap)
    cpdef int collapse(self) except -1
    cdef void _renderAppend(self, str s) except *
    cpdef str render(self, bint lengths=*, lengthFormat=*, Taxa.Map taxaMap=*,
      file outFile=*)

cdef class Node:
    cdef readonly Tree tree
    cdef readonly Ring ring
    cdef Taxon _taxon
    cdef int _degree
    cdef public object aux

    cdef Taxon getTaxon(self)
    cdef void setTaxon(self, Taxon taxon) except *
    # property taxon
    cpdef int _degreeGet(self, bint calculate=*) except -1
    cdef int getDegree(self)
    # property degree
    cdef rrender(self, Edge via, bint lengths, lengthFormat, Taxa.Map taxaMap,
      bint zeroLength, bint noLength)
    cpdef int separation(self, Node other)

cdef class Edge:
    cdef readonly Tree tree
    cdef readonly Ring ring
    cdef public double length
    cdef public object aux

    cpdef attach(self, Node nodeA, Node nodeB)
    cpdef detach(self)

cdef class Ring:
    cdef readonly Node node
    cdef readonly Edge edge
    cdef readonly Ring other
    cdef readonly Ring next
    cdef readonly Ring prev
    cdef public object aux

    cpdef siblings(self) # Iterator.

    cdef Node _minTaxon(self, Taxa.Map taxaMap)
    cdef Node _someLeaf(self)
    cdef Node _canonize(self, Taxa.Map taxaMap)
    cdef void _collapsable(self, list collapsable, list clampable) except *
    cdef void _collapse(self)
    cdef int _separation(self, Node other, int sep)
