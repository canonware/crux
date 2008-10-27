# Forward declarations.
cdef class Tree
cdef class Node
cdef class Edge
cdef class Ring

from CTMatrix cimport CTMatrix
cimport Newick
from TaxonMap cimport TaxonMap

import __builtin__

import random
import re

class _NewickParser(Newick.Parser):
    def __init__(self, tree, taxonMap, newickAutoMap=False):
        self._tree = tree
        self._taxonMap = taxonMap
        self._taxonStack = []
        self._newickAutoMap = newickAutoMap

    # Overridden method.
    def parse(self, str input):
        if not Newick.Parser.parse(self, input):
            if len(self._taxonStack) > 0:
                self._tree.baseSet(self._taxonStack[0])
            rVal = False
        else:
            rVal = True

        return rVal

    # Overridden method.
    def openParenAccept(self):
        self._taxonStack.insert(0, None)

    # Overridden method.
    def closeParenAccept(self):
        # Create an internal node, and join the top nodes to it.
        cnt = 0
        for elm in self._taxonStack:
            if elm == None:
                break
            elif type(elm) == Node:
                cnt += 1
        if cnt < 2:
            # Not enough neighbors on stack to join nodes together.  Remove open
            # paren (represented as None) from stack.
            self._taxonStack.remove(None)
        else:
            # Create new node.
            nnode = Node(self._tree)

            # Iteratively connect neighboring nodes to nnode.
            for i in xrange(cnt):
                if type(self._taxonStack[0]) == float:
                    length = self._taxonStack.pop(0)
                else:
                    length = 0.0

                n = self._taxonStack.pop(0)
                e = Edge(self._tree)
                e.attach(nnode, n)
                e.lengthSet(length)

            # Pop paren (None).
            self._taxonStack.pop(0)

            # Push nnode onto the stack.
            self._taxonStack.insert(0, nnode)

    # Helper method, called by {root,leaf}LabelAccept().
    def _labelAccept(self):
        if self._taxonMap.indGet(self.token()) != None:
            # Taxon mapping defined.
            val = self._taxonMap.indGet(self.token())
        else:
            # No taxon mapping defined; try to convert the label to an
            # integer.
            try:
                val = int(self.token())
                self._taxonMap.map(self.token(), val)
            except __builtin__.ValueError:
                if self._newickAutoMap:
                    # Create a new mapping.
                    val = self._taxonMap.ntaxaGet()
                    self._taxonMap.map(self.token(), val)
                else:
                    # Failed conversion.
                    raise Tree.ValueError, \
                          "At offset %d: No mapping for '%s'" \
                          % (self.offset(), self.token())

        # Create a new node and push it onto the stack.
        nnode = Node(self._tree)
        nnode.taxonNumSet(val)
        self._taxonStack.insert(0, nnode)

    # Overridden method.
    def rootLabelAccept(self):
        if len(self.token()) > 0:
            self._labelAccept()

    # Overridden method.
    def leafLabelAccept(self):
        self._labelAccept()

    # Overridden method.
    def lengthAccept(self):
        self._taxonStack.insert(0, float(self.token()))

    # Overridden method.
    def semicolonAccept(self):
        # If there is an internal node with only two neighbors on the stack,
        # splice it out of the tree.
        if len(self._taxonStack) != 0:
            if type(self._taxonStack[0]) == float:
                length = self._taxonStack.pop(0)
            else:
                length = None

            n = self._taxonStack[0]
            if n.degree() == 2:
                self._taxonStack.pop(0)
                # Get rings.
                ringA = n.ring()
                ringB = ringA.next()
                # Get neighboring nodes.
                nodeA = ringA.other().node()
                nodeB = ringB.other().node()
                # Detach edges.
                ringA.edge().detach()
                ringB.edge().detach()
                # Attach neighbors.
                e = Edge(self._tree)
                e.attach(nodeA, nodeB)
                if length != None:
                    e.lengthSet(length)
                else:
                    e.lengthSet(ringA.edge().lengthGet()
                                + ringB.edge().lengthGet())

                # Push a node back onto the stack.
                self._taxonStack.insert(0, nodeA)

# Default branch length callback function for random tree construction.
def _defaultRandomBranchCallback():
    return 1.0

# XXX Add some sort of taxon re-numbering method?
cdef class Tree:
    def __init__(self, with_=None, TaxonMap taxonMap=None,
      bint newickAutoMap=False, randomBranchCallback=None):
        self._sn = 0
        self._cacheSn = -1
        if type(with_) == int:
            self._randomNew(with_, taxonMap, randomBranchCallback)
        elif type(with_) == str or type(with_) == file:
            if taxonMap == None:
                taxonMap = TaxonMap()
            self._taxonMap = taxonMap
            self._newickNew(with_, newickAutoMap)
        else:
            if taxonMap == None:
                taxonMap = TaxonMap()
            self._taxonMap = taxonMap

    def _randomNew(self, ntaxa, taxonMap, randomBranchCallback):
        if taxonMap == None:
            self._taxonMap = TaxonMap()
            for i in xrange(ntaxa):
                self._taxonMap.map("T%d" % i, i)
        else:
            self._taxonMap = taxonMap

        # By default, generate branches of length 1.0.
        if randomBranchCallback == None:
            randomBranchCallback = _defaultRandomBranchCallback

        if ntaxa == 0:
            return

        # Create the first taxon.
        nodeA = Node(self)
        nodeA.taxonNumSet(0)
        self.baseSet(nodeA)
        if ntaxa == 1:
            return
        # Create the initial two-taxon tree that other taxa will be added to.
        nodeB = Node(self)
        nodeB.taxonNumSet(1)
        edge = Edge(self)
        edges = [edge]
        edge.lengthSet(randomBranchCallback())
        edge.attach(nodeA, nodeB)

        # Use random sequential addition to attach the remaning taxa.
        for i in xrange(2, ntaxa):
            # Pick an edge to bisect and add this taxon to.
            edgeA = edges[random.randint(0, len(edges)-1)]

            # Attach a new taxon node to a new internal node.
            nodeA = Node(self)
            nodeB = Node(self)
            nodeB.taxonNumSet(i)
            edgeB = Edge(self)
            edges.append(edgeB)
            edgeB.lengthSet(randomBranchCallback())
            edgeB.attach(nodeA, nodeB)

            (ringA, ringB) = edgeA.rings()
            nodeB = ringA.node()
            nodeC = ringB.node()
            edgeA.detach()
            edgeA.attach(nodeA, nodeB)

            edgeB = Edge(self)
            edges.append(edgeB)
            edgeB.lengthSet(randomBranchCallback())
            edgeB.attach(nodeA, nodeC)

    def _newickNew(self, input, newickAutoMap):
        parser = _NewickParser(self, self._taxonMap, newickAutoMap)
        return parser.parse(input)

    def _dup(self, newTree, node, prevRing):
        newNode = Node(newTree)
        taxonNum = node.taxonNumGet()
        if taxonNum != None:
            newNode.taxonNumSet(taxonNum)

        i = 0
        degree = node.degree()
        ring = node.ring()
        while (i < degree):
            if ring != prevRing:
                otherRing = ring.other()
                newOtherNode = self._dup(newTree, otherRing.node(), otherRing)

                newEdge = Edge(newTree)
                newEdge.attach(newNode, newOtherNode)

            ring = ring.next()
            i += 1

        return newNode

    cpdef dup(self):
        newTree = Tree()
        newNode = self._dup(newTree, self.baseGet(), None)
        newTree.baseSet(newNode)

    # XXX Make a property.
    cpdef taxonMapGet(self):
        return self._taxonMap

    cpdef rf(self, Tree other):
        if type(other) == Tree:
            return self._rfPair(other)
        else:
            return self._rfSequence(other)

    cdef _recacheRecurse(self, Ring ring):
        cdef Ring r

        if ring._node._taxonNum != -1:
            self._cachedNtaxa += 1

        for r in ring.siblings():
            self._cachedNedges += 1
            self._recacheRecurse(r.other())

    cdef _recache(self):
        cdef Ring ring, r

        self._cachedNtaxa = 0
        self._cachedNedges = 0

        if self._base != None:
            ring = self._base._ring
            if ring != None:
                if ring._node._taxonNum != -1:
                    self._cachedNtaxa += 1
                for r in ring:
                    self._cachedNedges += 1
                    self._recacheRecurse(r.other())

        self._cacheSn = self._sn

    # XXX Make a property.
    cpdef ntaxaGet(self):
        if self._cacheSn != self._sn:
            self._recache()
        return self._cachedNtaxa

    # XXX Make a property.
    cpdef nedgesGet(self):
        if self._cacheSn != self._sn:
            self._recache()
        return self._cachedNedges

    # XXX Make a property.
    cpdef Node baseGet(self):
        return self._base
    cpdef baseSet(self, Node base):
        self._base = base
        self._sn += 1

    cpdef canonize(self):
        pass # XXX
        self._sn += 1

    cpdef collapse(self):
        pass # XXX
        self._sn += 1

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
        # Make sure that taxon maps are identical.
        if not self._taxonMap.equal(cTMatrix.taxonMapGet()):
            raise Tree.ValueError(
                "TaxonMaps for Tree and CTMatrix must be equal")

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

    cpdef render(self, bint rooted=False, bint labels=False, bint lengths=False,
      lengthFormat="%.5e", outFile=None):
        cdef ret
        cdef Node n, neighbor
        cdef Ring ring

        # Set up for sending output to either a string or a file.
        if outFile == None:
            callback = self._stringRenderCallback
            self._renderTarget = ""
        else:
            callback = self._fileRenderCallback
            self._renderTarget = outFile

        # Render.
        n = self._base
        if n != None:
            if n.taxonNumGet() != None:
                # Leaf node.  If this node's neighbor is an internal node, start
                # rendering with it, in order to unroot the tree.
                ring = n.ring()
                if ring != None:
                    neighbor = ring.other().node()
                    if neighbor.taxonNumGet() == None:
                        if rooted:
                            n.rrender(None, self._taxonMap, labels, lengths,
                              lengthFormat, callback)
                        else:
                            # Start with the internal node.
                            neighbor.rrender(None, self._taxonMap, labels,
                              lengths, lengthFormat, callback)
                        callback(";")
                    else:
                        # This tree only has two taxa; start with the tree base.
                        callback("(")
                        n.rrender(None, self._taxonMap, labels, lengths,
                                  lengthFormat, callback, twoTaxa=True)
                        callback(");")
                else:
                    # There is only one node in the tree.
                    n.rrender(None, self._taxonMap, labels, lengths,
                              lengthFormat, callback)
                    callback(";")
            else:
                # Internal node.
                n.rrender(None, self._taxonMap, labels, lengths, lengthFormat,
                          callback)
                callback(";")
        else:
            callback(";")

        # Clean up and set ret according to where the output was sent.
        if outFile == None:
            ret = self._renderTarget
            self._renderTarget = None
        else:
            callback("\n")
            ret = None
            self._renderTarget = None
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
        self._taxonNum = -1

    # XXX Make a property.
    cpdef Tree tree(self):
        return self._tree

    # XXX Make a property.
    cpdef int taxonNumGet(self):
        return self._taxonNum
    cpdef taxonNumSet(self, int taxonNum):
        self._taxonNum = taxonNum
        self._tree._sn += 1

    # XXX Make a property.
    cpdef Ring ring(self):
        return self._ring

    # XXX Make a property.
    cpdef int degree(self):
        cdef int ret
        cdef Ring ring

        if self._ring == None:
            return 0
        ret = 1
        ring = self._ring._next
        while ring != self._ring:
            ret += 1
            ring = ring._next
        return ret

    cpdef rrender(self, Node prev, TaxonMap taxonMap, bint labels,
      bint lengths, lengthFormat, callback, bint twoTaxa=False):
        cdef bint did_something = False
        cdef bint did_paren = False
        cdef object taxonLabel, m
        cdef Ring ring
        cdef Node neighbor

        if self._taxonNum != -1:
            # Leaf node.
            if labels:
                # Protect special characters, if necessary.
                taxonLabel = taxonMap.labelGet(self._taxonNum)
                m = re.compile("[^ ()[\]':;,]*[ ()[\]':;,]").match(taxonLabel)
                if m:
                    taxonLabel = re.compile("'").sub("''", taxonLabel)
                    callback("'%s'" % taxonLabel)
                else:
                    callback("%s" % taxonLabel)
            else:
                callback("%d" % self.taxonNumGet())

            if lengths:
                ring = self._ring
                if ring != None:
                    if twoTaxa:
                        # This tree only has two taxa; take care not to double
                        # the branch length.
                        callback((":" + lengthFormat) \
                                 % (ring._edge._length / 2))
                    else:
                        callback((":" + lengthFormat) \
                                 % (ring._edge._length))
            did_something = True

        # Iterate through neighbors.
        ring = self._ring
        if ring != None:
            for i in xrange(self.degree()):
                # Get the node on the other end of the edge.  If it isn't prev,
                # recurse.
                neighbor = ring.other()._node
                if neighbor != prev:
                    if did_something:
                        callback(",")
                    elif not did_paren:
                        callback("(")
                        did_paren = True
                        did_something = True

                    neighbor.rrender(self, taxonMap, labels, lengths,
                                     lengthFormat, callback, twoTaxa)

                    if lengths:
                        if neighbor._taxonNum == -1:
                            callback((":" + lengthFormat) \
                                     % ring._edge._length)

                ring = ring._next

            if did_paren:
                callback(")")

    """Compute the number of edges that separate self and other."""
    cpdef int separation(self, Node other):
        cdef int ret
        cdef Ring ring, r

        if other is self:
            return 0

        ring = self._ring
        if ring != None:
            for r in ring:
                ret = r.other()._separation(other, 1)
                if ret != -1:
                    return ret
        return -1

# XXX Make methods into properties.
cdef class Edge:
    def __init__(self, Tree tree):
        self._tree = tree
        self._length = 0.0
        self._ringA = Ring(self)
        self._ringB = Ring(self)

    cpdef Tree tree(self):
        return self._tree

    # XXX Directly expose both rings as attributes to avoid tuple creation.
    cpdef rings(self):
        return (self._ringA, self._ringB)

    cpdef float lengthGet(self):
        return self._length

    cpdef lengthSet(self, float length):
        self._length = length
        self._tree._sn += 1

    cpdef attach(self, Node nodeA, Node nodeB):
        cdef Ring ring, nRing, pRing

        assert self._ringA._node == None
        assert self._ringB._node == None
        assert nodeA.separation(nodeB) == -1

        ring = self._ringA
        ring._node = nodeA
        nRing = nodeA._ring
        if nRing != None:
            pRing = nRing._prev
            ring._next = nRing
            ring._prev = pRing
            nRing._prev = ring
            pRing._next = ring
        nodeA._ring = ring

        ring = self._ringB
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
        assert nodeA.separation(nodeB) == 1

    cpdef detach(self):
        cdef Ring ring, nRing, pRing
        cdef Node node

        assert type(self._ringA._node) == Node
        assert type(self._ringB._node) == Node
        assert self._ringA._node.separation(self._ringB._node) == 1

        for ring in (self._ringA, self._ringB):
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

# XXX Make methods into properties.
cdef class Ring:
    def __init__(self, Edge edge):
        self._node = None
        self._edge = edge
        self._next = self
        self._prev = self

    """Iterate over "sibling" ring, including self."""
    def __iter__(self):
        return _RingIterHelper(self, True)

    """Iterate over "sibling" ring, excluding self."""
    cpdef siblings(self):
        return _RingIterHelper(self, False)

    cpdef Tree tree(self):
        return self._edge._tree

    cpdef Node node(self):
        return self._node

    cpdef Edge edge(self):
        return self._edge

    cpdef Ring other(self):
        cdef Ring ret = self._edge._ringA
        if ret != self:
            return ret
        return self._edge._ringB

    cpdef Ring next(self):
        return self._next

    cpdef Ring prev(self):
        return self._prev

    cdef int _separation(self, Node other, int sep):
        cdef int ret
        cdef Ring r

        if self._node is other:
            return sep

        for r in self.siblings():
            ret = r.other()._separation(other, sep + 1)
            if ret != -1:
                return ret

        return -1
