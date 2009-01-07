import Crux.Exception

class Exception(Crux.Exception.Exception):
    pass

import exceptions

class Malformed(Exception, exceptions.SyntaxError):
    def __init__(self, str):
        self._str = str

    def __str__(self):
        return self._str

import random
import re
import sys
import weakref

from Crux.CTMatrix cimport CTMatrix
cimport Crux.Newick as Newick
from Crux.Taxa cimport Taxon
cimport Crux.Taxa as Taxa

import Crux.Config

global __name__

# Forward declarations.
cdef class Tree
cdef class Node
cdef class Edge
cdef class Ring

cdef class Tree:
    def __init__(self, with_=None, Taxa.Map taxaMap=None, bint rooted=True):
        self._taxa = weakref.WeakValueDictionary()
        self._nodes = weakref.WeakKeyDictionary()
        self._edges = weakref.WeakKeyDictionary()
        self._sn = 0
        self._cacheSn = -1
        self.rooted = rooted
        if type(with_) == int:
            self._randomNew(with_, taxaMap)
        elif type(with_) == str:
            self._newickNew(with_, taxaMap)
        else:
            assert with_ is None

    cdef void _randomNew(self, int ntaxa, Taxa.Map taxaMap) except *:
        cdef Node nodeA, nodeB, nodeC
        cdef Edge edge, edgeA, edgeB
        cdef Ring ringA, ringB
        cdef list edges
        cdef int i

        if taxaMap is None:
            taxaMap = Taxa.Map()
            for 0 <= i < ntaxa:
                taxaMap.map(Taxa.get("T%d" % i), i)

        # Create the root.
        nodeA = Node(self)
        self.base = nodeA
        if ntaxa == 0:
            return

        # Attach the first taxon.
        nodeB = Node(self)
        nodeB.taxon = taxaMap.taxonGet(0)
        edge = Edge(self)
        edges = [edge]
        edge.attach(nodeA, nodeB)

        # Use random sequential addition to attach the remaning taxa.
        for 1 <= i < ntaxa:
            # Pick an edge to bisect and add this taxon to.
            edgeA = <Edge>(edges[random.randint(0, (i - 1) * 2)])

            # Attach a new taxon node to a new internal node.
            nodeA = Node(self)
            nodeB = Node(self)
            nodeB.taxon = taxaMap.taxonGet(i)
            edgeB = Edge(self)
            edges.append(edgeB)
            edgeB.attach(nodeA, nodeB)

            ringA = edgeA._ring
            ringB = ringA._other
            nodeB = ringA._node
            nodeC = ringB._node
            edgeA.detach()
            edgeA.attach(nodeA, nodeB)

            edgeB = Edge(self)
            edges.append(edgeB)
            edgeB.attach(nodeA, nodeC)

        # Convert to unrooted if necessary.
        if not self.rooted:
            self.rooted = True
            self.deroot()

    cdef void _newickNew(self, str input, Taxa.Map taxaMap) except *:
        cdef Newick.Parser parser

        parser = Newick.Parser(self, taxaMap)
        parser.parse(input)

        # Convert to unrooted if necessary.
        if not self.rooted:
            self.rooted = True
            self.deroot()

    cdef Node _dup(self, Tree newTree, Node node, Ring prevRing):
        cdef Node newNode, newOtherNode
        cdef Taxon taxon
        cdef int i, degree
        cdef Ring ring, otherRing
        cdef Edge newEdge

        newNode = Node(newTree)
        taxon = node._taxon
        if taxon is not None:
            newNode.taxon = taxon

        i = 0
        degree = node._degreeGet()
        ring = node._ring
        while i < degree:
            if ring is not prevRing:
                otherRing = ring._other
                newOtherNode = self._dup(newTree, otherRing._node, otherRing)

                newEdge = Edge(newTree)
                newEdge.attach(newNode, newOtherNode)

            ring = ring._next
            i += 1

        return newNode

    cpdef Tree dup(self):
        cdef Tree newTree
        cdef Node newBase

        newTree = Tree()
        newBase = self._dup(newTree, self._base, None)
        newTree.base = newBase
        newTree.rooted = self.rooted

        return newTree

    property taxa:
        """
            Alphabetized list of taxa.  If taxa are removed, this property may
            not be immediately updated.
        """
        def __get__(self):
            cdef list ret
            ret = self._taxa.keys()
            ret.sort()
            return ret

    property nodes:
        """
            List of nodes.  If nodes are removed, this property may not be
            immediately updated.
        """
        def __get__(self):
            cdef list ret
            ret = self._nodes.keys()
            return ret

    property edges:
        """
            List of edges.  If edges are removed, this property may not be
            immediately updated.
        """
        def __get__(self):
            cdef list ret
            ret = self._edges.keys()
            return ret

    cpdef rf(self, Tree other):
        if type(other) == Tree:
            return self._rfPair(other)
        else:
            return self._rfSequence(other)

    cdef void _recacheRecurse(self, Ring ring):
        cdef Ring r
        cdef Node node

        node = ring._node
        if node._taxon is not None:
            self._cachedNtaxa += 1
        self._cachedNnodes += 1

        node._degree = 1
        for r in ring.siblings():
            self._cachedNedges += 1
            self._recacheRecurse(r._other)
            node._degree += 1

    cdef void _recache(self):
        cdef Ring ring, r
        cdef Node node

        self._cachedNtaxa = 0
        self._cachedNnodes = 0
        self._cachedNedges = 0

        if self._base != None:
            node = self._base
            if node._taxon is not None:
                self._cachedNtaxa += 1
            self._cachedNnodes += 1
            node._degree = 0
            ring = node._ring
            if ring != None:
                for r in ring:
                    self._cachedNedges += 1
                    self._recacheRecurse(r._other)
                    node._degree += 1

        self._cacheSn = self._sn

    property ntaxa:
        def __get__(self):
            if self._cacheSn != self._sn:
                self._recache()
            return self._cachedNtaxa

    property nnodes:
        def __get__(self):
            if self._cacheSn != self._sn:
                self._recache()
            return self._cachedNnodes

    property nedges:
        def __get__(self):
            if self._cacheSn != self._sn:
                self._recache()
            return self._cachedNedges

    property base:
        def __get__(self):
            return self._base
        def __set__(self, Node base):
            self._base = base
            self._sn += 1

    cpdef deroot(self):
        cdef Node node
        cdef Edge edge
        cdef Ring ring
        cdef float removedLength

        if not self.rooted:
            return
        self.rooted = False

        ring = self._base._ring
        if ring is not None:
            ring = ring._other
            node = ring._node
            edge = ring._edge
            edge.detach()
            if node._degreeGet(calculate=True) == 2:
                # Detaching the root left an internal node with two edges.
                # Splice the node out.
                ring = node._ring
                self.base = ring._other._node
                edge = ring._edge
                removedLength = edge._length
                edge.detach()
                ring = node._ring
                edge = ring._edge
                node = ring._other._node
                edge.detach()
                edge.attach(self._base, node)
                edge.length = edge._length + removedLength
                # Change ring header/order in order to avoid disturbing
                # canonical order.
                self._base._ring = self._base._ring._next
            else:
                self.base = node
        else:
            # Only the root node exists, so discard the whole tree.
            self.base = None

    cpdef canonize(self, Taxa.Map taxaMap):
        cdef Node base, node
        cdef Ring ring

        if self.rooted:
            base = self._base
            assert base is not None
            assert base._degreeGet(calculate=True) <= 1
        else:
            # Find the minimum taxon before canonizing, and set it as the tree
            # base.  This is critical to correct results, since starting from
            # any other location in the tree will cause per-subtree minimum
            # taxa to be determined from the wrong perspective.
            node = self._base
            if node is None:
                return
            ring = node._ring
            if ring is None:
                return
            base = ring._minTaxon(taxaMap)
            node = ring._other._minTaxon(taxaMap)
            if base._taxon is None or \
              taxaMap.indGet(node._taxon) < taxaMap.indGet(base._taxon):
                base = node
            self.base = base

        ring = base._ring
        if ring is None:
            return
        ring._other._canonize(taxaMap)

    cpdef int collapse(self) except -1:
        cdef list collapsable, clampable
        cdef Node minTaxon
        cdef Ring ring, r
        cdef Edge edge

        if self._base is None:
            return 0

        # If the base is an internal node, move the base in order to keep from
        # losing the tree, should the current base node be removed from the
        # tree.
        if self._base._degreeGet(calculate=True) > 1:
            self.base = self._base._ring._other._someLeaf()
        assert self._base._degreeGet() <= 1

        # Generate a list of collapsable edges (but actually store rings in the
        # list, in order to be able to tell which end of the edge is closer to
        # the tree base).
        #
        # Also generate a list of clampablel edge.  Leaf edges cannot be
        # removed, but their lengths can be clamped at 0.0.
        ring = self._base._ring
        if ring is None:
            return 0
        collapsable = []
        clampable = []
        for r in ring:
            r._other._collapsable(collapsable, clampable)

        # Collapse edges.
        for r in collapsable:
            r._collapse()

        # Clamp leaf edge lengths.
        for edge in clampable:
            edge.length = 0.0

        return len(collapsable)

    cpdef tbr(self, Edge bisect, Edge reconnectA, Edge reconnectB):
        pass # XXX
        self._sn += 1

    # XXX Make a property.
    cpdef int tbrNNeigbhorsGet(self):
        pass # XXX

    # XXX Make a property.
    cpdef tbrNeighborGet(self, int neighbor):
        pass # XXX

    cpdef nni(self, Edge edge, Edge reconnectA, Edge reconnectB):
        pass # XXX
        self._sn += 1

    # XXX Make a property.
    cpdef nniNNeigbhorsGet(self):
        pass # XXX

    # XXX Make a property.
    cpdef nniNeighborGet(self, int neighbor):
        pass # XXX

    # XXX Rename mp* methods to reflect that these implement *Fitch* parsimony
    # *scoring*.  *Maximum* parsimony is a different concept.
    cpdef mpPrepare(self, CTMatrix cTMatrix, bint elimUninformative=True):
        cdef Taxon taxon

        # Make sure that cTMatrix.taxaMap is compatible.
        if cTMatrix.taxaMap.ntaxa != len(self._taxa):
            raise Tree.ValueError(
                "Taxa.Map for Tree and CTMatrix must be equal")
        for taxon in self._taxa:
            if cTMatrix.taxaMap.indGet(taxon) == -1:
                raise Tree.ValueError(
                  "Taxa.Map for CTMatrix does not contain taxon: %s" %
                  taxon.label)

        self._mpPrepare(cTMatrix, elimUninformative)

    cpdef mpFinish(self):
        pass # XXX

    cpdef mp(self):
        pass # XXX

    # XXX Implement more sophisticated tree holding, such that TBR neighbors
    # can be merged into a general pool of held trees.
    cpdef tbrBestNeighbhorsMp(self, int maxHold=-1):
        pass # XXX

    cpdef tbrBetterNeighborsMp(self, int maxHold=-1):
        pass # XXX

    cpdef tbrAllNeighborsMp(self, int maxHold=-1):
        pass # XXX

    # XXX Make a property.
    cpdef nHeldGet(self):
        pass # XXX

    # XXX Make a property.
    cpdef heldGet(self, int i):
        pass # XXX

    cdef void _renderAppend(self, str s) except *:
        if self._renderList is not None:
            self._renderList.append(s)
        else:
            self._renderFile.write(s)

    cpdef str render(self, bint lengths=False, lengthFormat="%.7e",
      Taxa.Map taxaMap=None, file outFile=None):
        cdef str ret
        cdef Node n, neighbor
        cdef Ring ring
        cdef int degree

        if outFile is None:
            self._renderList = []
            self._renderFile = None
        else:
            self._renderList = None
            self._renderFile = outFile

        # Render.
        n = self._base
        if n != None:
            if self.rooted:
#                self._renderAppend("[&r] ")
                if n._taxon is not None:
                    raise Malformed("Root is labeled")

                degree = n._degreeGet()
                if degree > 1:
                    raise Malformed("Root is an internal node")

                if degree == 1:
                    ring = n._ring
                    assert ring != None
                    neighbor = ring._other._node
                    neighbor.rrender(n, lengths, lengthFormat, taxaMap)
            else: # Unrooted tree.
#                self._renderAppend("[&u] ")
                degree = n._degreeGet()
                if degree == 0:
                    # There is only one node in the tree.
                    n.rrender(None, lengths, lengthFormat, taxaMap)
                elif degree == 1:
                    # Leaf node.  If this node's neighbor is an internal node,
                    # start rendering with it, in order to unroot the tree.
                    ring = n._ring
                    assert ring != None
                    neighbor = ring._other._node
                    if neighbor._degreeGet() > 1:
                        # Start with the internal node.
                        neighbor.rrender(None, lengths, lengthFormat, taxaMap,
                          noLength=True)
                    else:
                        # This tree only has two taxa; start with the tree
                        # base.
                        self._renderAppend("(")
                        n.rrender(neighbor, lengths, lengthFormat, taxaMap)
                        self._renderAppend(",")
                        neighbor.rrender(n, lengths, lengthFormat, taxaMap,
                          zeroLength=True)
                        self._renderAppend(")")
                else:
                    # Internal node.
                    n.rrender(None, lengths, lengthFormat, taxaMap,
                      noLength=True)

        self._renderAppend(";")

        if self._renderList is not None:
            ret = "".join(self._renderList)
            self._renderList = None
        else:
            self._renderFile.write("\n")
            ret = None
            self._renderFile = None
        return ret

    # Callback method that is used by the render method for recursive rendering
    # of the tree in Newick format.
    def _stringRenderCallback(self, string):
        # Append string to previous strings that were passed to this callback.
        self._renderTarget = "%s%s" % (self._renderTarget, string)

    # Callback method that is used by the render method for recursive rendering
    # of the tree in Newick format.
    def _fileRenderCallback(self, string):
        # Print string to self._renderTarget.
        self._renderTarget.write(string)

cdef class Node:
    def __init__(self, Tree tree):
        self._tree = tree
        self._ring = None
        self._taxon = None

        tree._nodes[self] = None

    property tree:
        def __get__(self):
            return self._tree

    property taxon:
        def __get__(self):
            return self._taxon
        def __set__(self, Taxon taxon):
            if __debug__:
                if taxon is not self._taxon and taxon in self._tree._taxa \
                  and self._tree._base is not None:
                    node = self._tree._taxa[taxon]
                    if node.separation(self._tree._base) != -1:
                        raise Malformed("Taxon already in use: %r" % \
                          taxon.label)

            if self._taxon is not None:
                self._tree._taxa.pop(self._taxon)
            self._tree._taxa[taxon] = self
            self._taxon = taxon

    property ring:
        def __get__(self):
            return self._ring

    cpdef int _degreeGet(self, bint calculate=False) except -1:
        """
            Return the number of attached edges.  If the node is reachable via
            the tree base (self.separation(self.tree.base) != -1), it is
            possible to use a cached value.  For an unreachable node, or in
            order to avoid computing the cache, set calculate to True.
        """
        cdef int ret
        cdef Ring ring

        if not calculate:
            assert self._tree._base != None
#            assert self.separation(self._tree._base) != -1 # Expensive.

            if self._tree._cacheSn != self._tree._sn:
                self._tree._recache()

            return self._degree
        else:
            ret = 0
            if self._ring != None:
                for ring in self._ring:
                    ret += 1
            return ret

    property degree:
        def __get__(self):
            return self._degreeGet(calculate=False)

    cpdef rrender(self, Node prev, bint lengths, lengthFormat, Taxa.Map taxaMap,
      bint zeroLength=False, bint noLength=False):
        cdef bint did_paren = False
        cdef str label
        cdef object m
        cdef Ring ring, r
        cdef Node neighbor
        cdef int degree

        # Iterate through neighbors.
        ring = self._ring
        if ring != None:
            for r in ring:
                # Get the node on the other end of the edge.  If it isn't prev,
                # recurse.
                neighbor = r._other._node
                if neighbor != prev:
                    if did_paren:
                        self._tree._renderAppend(",")
                    elif not did_paren:
                        self._tree._renderAppend("(")
                        did_paren = True

                    neighbor.rrender(self, lengths, lengthFormat, taxaMap)

            if did_paren:
                self._tree._renderAppend(")")

        # Render label.
        if self._taxon is not None:
            if taxaMap is not None:
                self._tree._renderAppend("%d" %
                  taxaMap.indGet(self._taxon))
            else:
                # Protect special characters, if necessary.
                label = self._taxon.label
                m = re.compile(r"[_()[\]':;,]").search(label)
                if m:
                    label = re.compile("'").sub("''", label)
                    self._tree._renderAppend("'%s'" % label)
                else:
                    if label.find(" ") != -1:
                        label = re.compile(" ").sub("_", label)
                    self._tree._renderAppend("%s" % label)

        # Render branch length.
        degree = self._degreeGet()
        if lengths and degree > 0:
            ring = self._ring
            assert ring != None
            if zeroLength:
                # This tree only has two taxa; take care not to double the
                # branch length.
                self._tree._renderAppend((":" + lengthFormat) % 0.0)
            elif not noLength:
                self._tree._renderAppend((":" + lengthFormat) % \
                  (ring._edge._length))

    cpdef int separation(self, Node other):
        """
            Compute the number of edges that separate self and other.
        """
        cdef int ret
        cdef Ring ring, r

        if other is self:
            return 0

        ring = self._ring
        if ring != None:
            for r in ring:
                ret = r._other._separation(other, 1)
                if ret != -1:
                    return ret
        return -1

cdef class Edge:
    def __init__(self, Tree tree):
        cdef Ring other

        self._tree = tree
        self._length = 0.0
        self._ring = Ring(self, None)
        other = Ring(self, self._ring)
        self._ring._other = other

        tree._edges[self] = None

    property tree:
        def __get__(self):
            return self._tree

    property ring:
        def __get__(self):
            return self._ring

    property length:
        def __get__(self):
            return self._length
        def __set__(self, float length):
            self._length = length
            self._tree._sn += 1

    cpdef attach(self, Node nodeA, Node nodeB):
        cdef Ring ring, nRing, pRing

        assert self._ring._node == None
        assert self._ring._other._node == None
#        assert nodeA.separation(nodeB) == -1 # Expensive.

        ring = self._ring
        ring._node = nodeA
        nRing = nodeA._ring
        if nRing != None:
            pRing = nRing._prev
            ring._next = nRing
            ring._prev = pRing
            nRing._prev = ring
            pRing._next = ring
        nodeA._ring = ring

        ring = self._ring._other
        ring._node = nodeB
        nRing = nodeB._ring
        if nRing != None:
            pRing = nRing._prev
            ring._next = nRing
            ring._prev = pRing
            nRing._prev = ring
            pRing._next = ring
        nodeB._ring = ring

        self._tree._sn += 1
#        assert nodeA.separation(nodeB) == 1 # Expensive.

    cpdef detach(self):
        cdef Ring ring, nRing, pRing
        cdef Node node

        assert type(self._ring._node) == Node
        assert type(self._ring._other._node) == Node
#        assert self._ring._node.separation(self._ring._other._node) == 1

        for ring in (self._ring, self._ring._other):
            node = ring._node
            nRing = ring._next
            if nRing == ring:
                node._ring = None
            else:
                if node._ring == ring:
                    node._ring = nRing
                pRing = ring._prev
                nRing._prev = pRing
                pRing._next = nRing
                ring._next = ring
                ring._prev = ring
            ring._node = None

        self._tree._sn += 1

cdef class _RingIterHelper:
    cdef Ring _start, _next

    def __init__(self, Ring ring, bint all):
        self._start = ring
        if all:
            self._next = None
        else:
            self._next = ring._next

    def __iter__(self):
        return self

    def __next__(self):
        cdef Ring ret

        if self._next == None:
            ret = self._start
            self._next = ret._next
        else:
            ret = self._next
            if ret == self._start:
                raise StopIteration
            self._next = ret._next
        return ret

cdef class Ring:
    def __init__(self, Edge edge, Ring other):
        self._node = None
        self._edge = edge
        self._other = other
        self._next = self
        self._prev = self

    def __iter__(self):
        """
            Iterate over "sibling" ring, including self.
        """
        return _RingIterHelper(self, True)

    cpdef siblings(self):
        """
            Iterate over "sibling" ring, excluding self.
        """
        return _RingIterHelper(self, False)

    cdef Node _minTaxon(self, Taxa.Map taxaMap):
        cdef Node ret, minTaxon
        cdef Ring r

        ret = self._node
        for r in self.siblings():
            minTaxon = r._other._minTaxon(taxaMap)
            if ret._taxon is None or \
              taxaMap.indGet(minTaxon._taxon) < taxaMap.indGet(ret._taxon):
                ret = minTaxon

        return ret

    cdef Node _someLeaf(self):
        cdef Ring next

        next = self._next
        if next is self:
            return self._node
        else:
            return next._other._someLeaf()

    cdef Node _canonize(self, Taxa.Map taxaMap):
        cdef Node ret, minTaxon, node, nodeOther
        cdef int degree
        cdef list rings
        cdef Ring r
        cdef Edge edge

        node = self._node
        ret = node

        degree = node._degreeGet(calculate=True)
        if degree > 1:
            rings = []

            # Iteratively canonize subtrees, keeping track of the minimum
            # taxon seen overall, as well as for each subtree.
            for r in self.siblings():
                minTaxon = r._other._canonize(taxaMap)
                if ret._taxon is None or \
                  taxaMap.indGet(minTaxon._taxon) < taxaMap.indGet(ret._taxon):
                    ret = minTaxon
                rings.append((taxaMap.indGet(minTaxon._taxon), r))

            # Sort according to per-subtree minimum taxa.
            rings.sort()

            # Detach and re-attach all edges, in reverse order.  This code
            # assumes that attaching inserts the edge at the head of the ring.
            for len(rings) > i >= 0:
                r = <Ring>rings[i][1]
                nodeOther = r._other._node
                edge = r._edge
                edge.detach()
                edge.attach(node, nodeOther)

            # Set the beginning of the ring to self.  This makes it easier for
            # external code to traverse a tree in canonical order.
            #
            # Note that there is no need to increment the tree sequence number,
            # since the edge detach/attach operations already do so.
            node._ring = self

        return ret

    cdef void _collapsable(self, list collapsable, list clampable) except *:
        cdef Ring ring
        cdef Edge edge

        for ring in self.siblings():
            ring._other._collapsable(collapsable, clampable)

        edge = self._edge
        if edge._length <= 0.0:
            if self._node._degreeGet() > 1 and \
              self._other._node._degreeGet() > 1:
                collapsable.append(self)
            else:
                # Leaf node.  Clamp length.
                clampable.append(edge)

    cdef void _collapse(self):
        cdef Ring rOther, rTemp
        cdef Edge edge, eTemp
        cdef Node node, nOther, nTemp

        # Collapse the edge that self is a part of.  At the end of this method,
        # self, edge, rOther, and nOther will have been removed from the tree.
        # Following is a diagram of how variables are related just before
        # detaching an eTemp.
        #
        #                       nTemp
        #                         |
        #                         |
        #                       [ring]
        #                         |
        #                         |
        #                       eTemp
        #                         |
        #                         |
        #                       rTemp
        #                      /
        #                     /
        # node--self--edge--rOther--nOther
        #       ^^^^^^^^^^^^^^^^^^^^^^^^^^
        #                Remove
        node = self._node
        edge = self._edge
        rOther = self._other
        nOther = rOther._node
        assert nOther._degreeGet(calculate=True) > 1

        rTemp = rOther._next
        while rTemp != rOther:
            eTemp = rTemp._edge
            nTemp = rTemp._other._node
            eTemp.detach()
            eTemp.attach(node, nTemp)
            rTemp = rOther._next

        assert nOther._degreeGet(calculate=True) == 1
        edge.detach()

    cdef int _separation(self, Node other, int sep):
        cdef int ret
        cdef Ring r

        if self._node is other:
            return sep

        for r in self.siblings():
            ret = r._other._separation(other, sep + 1)
            if ret != -1:
                return ret

        return -1

    property tree:
        def __get__(self):
            return self._edge._tree

    property node:
        def __get__(self):
            return self._node

    property edge:
        def __get__(self):
            return self._edge

    property other:
        def __get__(self):
            return self._other

    property next:
        def __get__(self):
            return self._next

    property prev:
        def __get__(self):
            return self._prev
