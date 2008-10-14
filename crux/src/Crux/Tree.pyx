from CTMatrix cimport CTMatrix
import NewickParser
from TaxonMap cimport TaxonMap

import Crux

import __builtin__

import random
import re

class _NewickParser(NewickParser.NewickParser):
    def __init__(self, tree, taxonMap, newickAutoMap=False):
        self._tree = tree
        self._taxonMap = taxonMap
        self._taxonStack = []
        self._newickAutoMap = newickAutoMap

    # Overridden method.
    def parse(self, input):
        if not NewickParser.NewickParser.parse(self, input):
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
                    raise Crux.Tree.ValueError, \
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

    cdef dup(self):
        newTree = Crux.Tree.Tree()
        newNode = self._dup(newTree, self.baseGet(), None)
        newTree.baseSet(newNode)

    # XXX Make a property.
    cdef taxonMapGet(self):
        pass # XXX

    cdef rf(self, Tree other):
        if type(other) == Tree:
            return self._rfPair(other)
        else:
            return self._rfSequence(other)

    # XXX Make a property.
    cdef ntaxaGet(self):
        pass # XXX

    # XXX Make a property.
    cdef nedgesGet(self):
        pass # XXX

    # XXX Make a property.
    cdef Node baseGet(self):
        pass # XXX
    cdef baseSet(self, Node base):
        pass # XXX

    cdef canonize(self):
        pass # XXX

    cdef collapse(self):
        pass # XXX

    cdef tbr(self, Edge bisect, Edge reconnectA, Edge reconnectB):
        pass # XXX

    # XXX Make a property.
    cdef int tbrNNeigbhorsGet(self):
        pass # XXX

    # XXX Make a property.
    cdef tbrNeighborGet(self, int neighbor):
        pass # XXX

    cdef nni(self, Edge edge, Edge reconnectA, Edge reconnectB):
        pass # XXX

    # XXX Make a property.
    cdef nniNNeigbhorsGet(self):
        pass # XXX

    # XXX Make a property.
    cdef nniNeighborGet(self, int neighbor):
        pass # XXX

    # XXX Rename mp* methods to reflect that these implement *Fitch* parsimony
    # *scoring*.  *Maximum* parsimony is a different concept.
    cdef mpPrepare(self, CTMatrix cTMatrix, bint elimUninformative=True):
        # Make sure that taxon maps are identical.
        if not self._taxonMap.equal(cTMatrix.taxonMapGet()):
            raise Crux.Tree.ValueError(
                "TaxonMaps for Tree and CTMatrix must be equal")

        self._mpPrepare(cTMatrix, elimUninformative)

    cdef mpFinish(self):
        pass # XXX

    cdef mp(self):
        pass # XXX

    # XXX Implement more sophisticated tree holding, such that TBR neighbors
    # can be merged into a general pool of held trees.
    cdef tbrBestNeighbhorsMp(self, int maxHold=-1):
        pass # XXX

    cdef tbrBetterNeighborsMp(self, int maxHold=-1):
        pass # XXX

    cdef tbrAllNeighborsMp(self, int maxHold=-1):
        pass # XXX

    # XXX Make a property.
    cdef nHeldGet(self):
        pass # XXX

    # XXX Make a property.
    cdef heldGet(self, int i):
        pass # XXX

    cdef render(self, bint rooted=False, bint labels=False, bint lengths=False,
      lengthFormat="%.5e", outFile=None):
        # Set up for sending output to either a string or a file.
        if outFile == None:
            callback = self._stringRenderCallback
            self._renderString = ""
        else:
            callback = self._fileRenderCallback
            self._renderOutFile = outFile

        # Render.
        n = self.baseGet()
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

        # Clean up and set rVal according to where the output was sent.
        if outFile == None:
            rVal = self._renderString
            self._renderString = None
        else:
            callback("\n")
            rVal = None
            self._renderOutFile = None
        return rVal

    # Callback method that is used by the render method for recursive rendering
    # of the tree in Newick format.
    def _stringRenderCallback(self, string):
        # Append string to previous strings that were passed to this callback.
        self._renderString = "%s%s" % (self._renderString, string)

    # Callback method that is used by the render method for recursive rendering
    # of the tree in Newick format.
    def _fileRenderCallback(self, string):
        # Print string to self._renderOutFile.
        self._renderOutFile.write(string)

cdef class Node:
    def __init__(self, tree):
        pass

    # XXX Make a property.
    cdef Tree tree(self):
        pass # XXX

    # XXX Make a property.
    cdef int taxonNumGet(self):
        pass # XXX
    cdef void taxonNumSet(self, int taxonNum):
        pass # XXX

    # XXX Make a property.
    cdef Ring ring(self):
        pass # XXX

    # XXX Make a property.
    cdef int degree(self):
        pass # XXX

    def rrender(self, prev, taxonMap, labels, lengths, lengthFormat, callback,
                twoTaxa=False):
        did_something = False
        did_paren = False

        if self.taxonNumGet() != None:
            # Leaf node.
            if labels:
                # Protect special characters, if necessary.
                taxonLabel = taxonMap.labelGet(self.taxonNumGet())
                m = re.compile("[^ ()[\]':;,]*[ ()[\]':;,]").match(taxonLabel)
                if m:
                    taxonLabel = re.compile("'").sub("''", taxonLabel)
                    callback("'%s'" % taxonLabel)
                else:
                    callback("%s" % taxonLabel)
            else:
                callback("%d" % self.taxonNumGet())

            if lengths:
                ring = self.ring()
                if ring != None:
                    if twoTaxa:
                        # This tree only has two taxa; take care not to double
                        # the branch length.
                        callback((":" + lengthFormat) \
                                 % (ring.edge().lengthGet() / 2))
                    else:
                        callback((":" + lengthFormat) \
                                 % (ring.edge().lengthGet()))
            did_something = True

        # Iterate through neighbors.
        ring = self.ring()
        if ring != None:
            for i in xrange(self.degree()):
                # Get the node on the other end of the edge.  If it isn't prev,
                # recurse.
                neighbor = ring.other().node()
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
                        if neighbor.taxonNumGet() == None:
                            callback((":" + lengthFormat) \
                                     % ring.edge().lengthGet())

                ring = ring.next()

            if did_paren:
                callback(")")

# XXX Make methods into properties.
cdef class Edge:
    def __init__(self, Tree tree):
        self._tree = tree
        self._length = 0.0
        self._ringA = Ring(self)
        self._ringB = Ring(self)

    cdef Tree tree(self):
        return self._tree

    # XXX Directly expose both rings as attributes to avoid tuple creation.
    cdef rings(self):
        return (self._ringA, self._ringB)

    cdef float lengthGet(self):
        return self._length

    cdef void lengthSet(self, float length):
        self._length = length

    cdef void attach(self, Node nodeA, Node nodeB):
        cdef Ring ring, nRing, pRing
        assert self._ringA._node == None
        assert self._ringB._node == None

        ring = self._ringA
        ring._node = nodeA
        nRing = nodeA._ring
        if nRing != None:
            pRing = nRing._prev
            ring._next = pRing
            ring._prev = nRing
            nRing._prev = ring
            pRing._next = ring
        nodeA._ring = ring

        ring = self._ringB
        ring._node = nodeB
        nRing = nodeB._ring
        if ring != None:
            pRing = nRing._prev
            ring._next = pRing
            ring._prev = nRing
            nRing._prev = ring
            pRing._next = ring
        nodeB._ring = ring

    cdef detach(self):
        assert type(self._ringA._node) == Node
        assert type(self._ringB._node) == Node

        self._ringA._node = None
        self._ringB._node = None

        # XXX

# XXX Make methods into properties.
cdef class Ring:
    cdef Tree tree(self):
        return self._edge.tree()

    cdef Node node(self):
        return self._node

    cdef Edge edge(self):
        return self._edge

    cdef Ring other(self):
        ret = self._edge._ringA
        if ret != self:
            return ret
        return self._edge._ringB

    cdef Ring next(self):
        return self._next

    cdef Ring prev(self):
        return self._prev
