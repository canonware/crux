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

            ringA = edgeA.ring
            ringB = ringA.other
            nodeB = ringA.node
            nodeC = ringB.node
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
        ring = node.ring
        while i < degree:
            if ring is not prevRing:
                otherRing = ring.other
                newOtherNode = self._dup(newTree, otherRing.node, otherRing)

                newEdge = Edge(newTree)
                newEdge.attach(newNode, newOtherNode)

            ring = ring.next
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

        node = ring.node
        if node._taxon is not None:
            self._cachedNtaxa += 1
        self._cachedNnodes += 1

        node._degree = 1
        for r in ring.siblings():
            self._cachedNedges += 1
            self._recacheRecurse(r.other)
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
            ring = node.ring
            if ring != None:
                for r in ring:
                    self._cachedNedges += 1
                    self._recacheRecurse(r.other)
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

        ring = self._base.ring
        if ring is not None:
            ring = ring.other
            node = ring.node
            edge = ring.edge
            edge.detach()
            if node._degreeGet(calculate=True) == 2:
                # Detaching the root left an internal node with two edges.
                # Splice the node out.
                ring = node.ring
                self.base = ring.other.node
                edge = ring.edge
                removedLength = edge._length
                edge.detach()
                ring = node.ring
                edge = ring.edge
                node = ring.other.node
                edge.detach()
                edge.attach(self._base, node)
                edge.length = edge._length + removedLength
                # Change ring header/order in order to avoid disturbing
                # canonical order.
                self._base.ring = self._base.ring.next
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
            ring = node.ring
            if ring is None:
                return
            base = ring._minTaxon(taxaMap)
            node = ring.other._minTaxon(taxaMap)
            if base._taxon is None or \
              taxaMap.indGet(node._taxon) < taxaMap.indGet(base._taxon):
                base = node
            self.base = base

        ring = base.ring
        if ring is None:
            return
        ring.other._canonize(taxaMap)

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
            self.base = self._base.ring.other._someLeaf()
        assert self._base._degreeGet() <= 1

        # Generate a list of collapsable edges (but actually store rings in the
        # list, in order to be able to tell which end of the edge is closer to
        # the tree base).
        #
        # Also generate a list of clampablel edge.  Leaf edges cannot be
        # removed, but their lengths can be clamped at 0.0.
        ring = self._base.ring
        if ring is None:
            return 0
        collapsable = []
        clampable = []
        for r in ring:
            r.other._collapsable(collapsable, clampable)

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
                    ring = n.ring
                    assert ring != None
                    neighbor = ring.other.node
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
                    ring = n.ring
                    assert ring != None
                    neighbor = ring.other.node
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
        self.tree = tree
        self.ring = None
        self._taxon = None

        tree._nodes[self] = None

    property taxon:
        def __get__(self):
            return self._taxon
        def __set__(self, Taxon taxon):
            if __debug__:
                if taxon is not self._taxon and taxon in self.tree._taxa \
                  and self.tree._base is not None:
                    node = self.tree._taxa[taxon]
                    if node.separation(self.tree._base) != -1:
                        raise Malformed("Taxon already in use: %r" % \
                          taxon.label)

            if self._taxon is not None:
                self.tree._taxa.pop(self._taxon)
            self.tree._taxa[taxon] = self
            self._taxon = taxon

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
            assert self.tree._base != None
#            assert self.separation(self.tree._base) != -1 # Expensive.

            if self.tree._cacheSn != self.tree._sn:
                self.tree._recache()

            return self._degree
        else:
            ret = 0
            if self.ring != None:
                for ring in self.ring:
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
        ring = self.ring
        if ring != None:
            for r in ring:
                # Get the node on the other end of the edge.  If it isn't prev,
                # recurse.
                neighbor = r.other.node
                if neighbor != prev:
                    if did_paren:
                        self.tree._renderAppend(",")
                    elif not did_paren:
                        self.tree._renderAppend("(")
                        did_paren = True

                    neighbor.rrender(self, lengths, lengthFormat, taxaMap)

            if did_paren:
                self.tree._renderAppend(")")

        # Render label.
        if self._taxon is not None:
            if taxaMap is not None:
                self.tree._renderAppend("%d" %
                  taxaMap.indGet(self._taxon))
            else:
                # Protect special characters, if necessary.
                label = self._taxon.label
                m = re.compile(r"[_()[\]':;,]").search(label)
                if m:
                    label = re.compile("'").sub("''", label)
                    self.tree._renderAppend("'%s'" % label)
                else:
                    if label.find(" ") != -1:
                        label = re.compile(" ").sub("_", label)
                    self.tree._renderAppend("%s" % label)

        # Render branch length.
        degree = self._degreeGet()
        if lengths and degree > 0:
            ring = self.ring
            assert ring != None
            if zeroLength:
                # This tree only has two taxa; take care not to double the
                # branch length.
                self.tree._renderAppend((":" + lengthFormat) % 0.0)
            elif not noLength:
                self.tree._renderAppend((":" + lengthFormat) % \
                  (ring.edge._length))

    cpdef int separation(self, Node other):
        """
            Compute the number of edges that separate self and other.
        """
        cdef int ret
        cdef Ring ring, r

        if other is self:
            return 0

        ring = self.ring
        if ring != None:
            for r in ring:
                ret = r.other._separation(other, 1)
                if ret != -1:
                    return ret
        return -1

cdef class Edge:
    def __init__(self, Tree tree):
        cdef Ring other

        self.tree = tree
        self._length = 0.0
        self.ring = Ring(self, None)
        other = Ring(self, self.ring)
        self.ring.other = other

        tree._edges[self] = None

    property length:
        def __get__(self):
            return self._length
        def __set__(self, float length):
            self._length = length
            self.tree._sn += 1

    cpdef attach(self, Node nodeA, Node nodeB):
        cdef Ring ring, nRing, pRing

        assert self.ring.node == None
        assert self.ring.other.node == None
#        assert nodeA.separation(nodeB) == -1 # Expensive.

        ring = self.ring
        ring.node = nodeA
        nRing = nodeA.ring
        if nRing != None:
            pRing = nRing.prev
            ring.next = nRing
            ring.prev = pRing
            nRing.prev = ring
            pRing.next = ring
        nodeA.ring = ring

        ring = self.ring.other
        ring.node = nodeB
        nRing = nodeB.ring
        if nRing != None:
            pRing = nRing.prev
            ring.next = nRing
            ring.prev = pRing
            nRing.prev = ring
            pRing.next = ring
        nodeB.ring = ring

        self.tree._sn += 1
#        assert nodeA.separation(nodeB) == 1 # Expensive.

    cpdef detach(self):
        cdef Ring ring, nRing, pRing
        cdef Node node

        assert type(self.ring.node) == Node
        assert type(self.ring.other.node) == Node
#        assert self.ring.node.separation(self.ring.other.node) == 1

        for ring in (self.ring, self.ring.other):
            node = ring.node
            nRing = ring.next
            if nRing == ring:
                node.ring = None
            else:
                if node.ring == ring:
                    node.ring = nRing
                pRing = ring.prev
                nRing.prev = pRing
                pRing.next = nRing
                ring.next = ring
                ring.prev = ring
            ring.node = None

        self.tree._sn += 1

cdef class _RingIterHelper:
    cdef Ring _start, _next

    def __init__(self, Ring ring, bint all):
        self._start = ring
        if all:
            self._next = None
        else:
            self._next = ring.next

    def __iter__(self):
        return self

    def __next__(self):
        cdef Ring ret

        if self._next == None:
            ret = self._start
            self._next = ret.next
        else:
            ret = self._next
            if ret == self._start:
                raise StopIteration
            self._next = ret.next
        return ret

cdef class Ring:
    def __init__(self, Edge edge, Ring other):
        self.node = None
        self.edge = edge
        self.other = other
        self.next = self
        self.prev = self

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

        ret = self.node
        for r in self.siblings():
            minTaxon = r.other._minTaxon(taxaMap)
            if ret._taxon is None or \
              taxaMap.indGet(minTaxon._taxon) < taxaMap.indGet(ret._taxon):
                ret = minTaxon

        return ret

    cdef Node _someLeaf(self):
        cdef Ring next

        next = self.next
        if next is self:
            return self.node
        else:
            return next.other._someLeaf()

    cdef Node _canonize(self, Taxa.Map taxaMap):
        cdef Node ret, minTaxon, node, nodeOther
        cdef int degree
        cdef list rings
        cdef Ring r
        cdef Edge edge

        node = self.node
        ret = node

        degree = node._degreeGet(calculate=True)
        if degree > 1:
            rings = []

            # Iteratively canonize subtrees, keeping track of the minimum
            # taxon seen overall, as well as for each subtree.
            for r in self.siblings():
                minTaxon = r.other._canonize(taxaMap)
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
                nodeOther = r.other.node
                edge = r.edge
                edge.detach()
                edge.attach(node, nodeOther)

            # Set the beginning of the ring to self.  This makes it easier for
            # external code to traverse a tree in canonical order.
            #
            # Note that there is no need to increment the tree sequence number,
            # since the edge detach/attach operations already do so.
            node.ring = self

        return ret

    cdef void _collapsable(self, list collapsable, list clampable) except *:
        cdef Ring ring
        cdef Edge edge

        for ring in self.siblings():
            ring.other._collapsable(collapsable, clampable)

        edge = self.edge
        if edge._length <= 0.0:
            if self.node._degreeGet() > 1 and \
              self.other.node._degreeGet() > 1:
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
        node = self.node
        edge = self.edge
        rOther = self.other
        nOther = rOther.node
        assert nOther._degreeGet(calculate=True) > 1

        rTemp = rOther.next
        while rTemp != rOther:
            eTemp = rTemp.edge
            nTemp = rTemp.other.node
            eTemp.detach()
            eTemp.attach(node, nTemp)
            rTemp = rOther.next

        assert nOther._degreeGet(calculate=True) == 1
        edge.detach()

    cdef int _separation(self, Node other, int sep):
        cdef int ret
        cdef Ring r

        if self.node is other:
            return sep

        for r in self.siblings():
            ret = r.other._separation(other, sep + 1)
            if ret != -1:
                return ret

        return -1
