from CTMatrix cimport CTMatrix
cimport Newick
from TaxonMap cimport TaxonMap
cimport Parsing

import Crux.Config
import Crux.Exception

class Exception(Crux.Exception.Exception):
    pass

import exceptions

class Malformed(Exception, exceptions.SyntaxError):
    def __init__(self, str):
        self._str = str

    def __str__(self):
        return self._str

import __builtin__
import random
import re
import sys

global __name__

# Forward declarations.
cdef class _NewickTree(Newick.Tree)
cdef class _NewickDescendantList(Newick.DescendantList)
cdef class _NewickSubtreeList(Newick.SubtreeList)
cdef class _NewickSubtree(Newick.Subtree)
cdef class _NewickLabel(Newick.Label)
cdef class _NewickParser(Newick.Parser)
cdef class _NewickParser(Newick.Parser)
cdef class Tree
cdef class Node
cdef class Edge
cdef class Ring

#===============================================================================
# Begin Newick tree construction support.
#

cdef class _NewickTree(Newick.Tree):
    "%extend Tree"

    cdef Node root

    def __init__(self, Parsing.Lr parser):
        Newick.Tree.__init__(self, parser)

        self.root = Node(self.parser._tree)

    # (A,B)C:4.2;
    #
    #   .
    #   |4.2
    #   C
    #  / \
    # A   B
    cpdef reduceDRB(self, Newick.DescendantList DescendantList,
      Newick.Label Label, Newick.TokenColon colon,
      Newick.TokenBranchLength branchLength, Newick.TokenSemicolon semicolon):
        "%accept"
        cdef Edge edge

        self.parser._labelNode((<_NewickDescendantList>DescendantList).node,
          <_NewickLabel>Label)
        edge = Edge(self.parser._tree)
        edge.lengthSet(float(branchLength.raw))
        edge.attach(self.root, (<_NewickDescendantList>DescendantList).node)

    # (A,B)C;
    #
    #   .
    #   |
    #   C
    #  / \
    # A   B
    cpdef reduceDR(self, Newick.DescendantList DescendantList,
      Newick.Label Label, Newick.TokenSemicolon semicolon):
        "%accept"
        cdef Edge edge

        self.parser._labelNode((<_NewickDescendantList>DescendantList).node,
          <_NewickLabel>Label)
        edge = Edge(self.parser._tree)
        edge.attach(self.root, (<_NewickDescendantList>DescendantList).node)

    # (A,B):4.2;
    #
    #   .
    #   |4.2
    #   .
    #  / \
    # A   B
    cpdef reduceDB(self, Newick.DescendantList DescendantList,
      Newick.TokenColon colon, Newick.TokenBranchLength branchLength,
      Newick.TokenSemicolon semicolon):
        "%accept"
        cdef Edge edge

        edge = Edge(self.parser._tree)
        edge.lengthSet(float(branchLength.raw))
        edge.attach(self.root, (<_NewickDescendantList>DescendantList).node)

    # A:4.2;
    #
    #   .
    #   |4.2
    #   A
    cpdef reduceRB(self, Newick.Label Label, Newick.TokenColon colon,
      Newick.TokenBranchLength branchLength, Newick.TokenSemicolon semicolon):
        "%accept"
        cdef Node node
        cdef Edge edge

        node = Node(self.parser._tree)
        self.parser._labelNode(node, <_NewickLabel>Label)
        edge = Edge(self.parser._tree)
        edge.lengthSet(float(branchLength.raw))
        edge.attach(self.root, node)

    # (A,B);
    #
    #   .
    #   |
    #   .
    #  / \
    # A   B
    cpdef reduceD(self, Newick.DescendantList DescendantList,
      Newick.TokenSemicolon semicolon):
        "%accept"
        cdef Edge edge

        edge = Edge(self.parser._tree)
        edge.attach(self.root, (<_NewickDescendantList>DescendantList).node)

    # A;
    #
    #   .
    #   |
    #   A
    cpdef reduceR(self, Newick.Label Label, Newick.TokenSemicolon semicolon):
        "%accept"
        cdef Node node
        cdef Edge edge

        node = Node(self.parser._tree)
        self.parser._labelNode(node, <_NewickLabel>Label)
        edge = Edge(self.parser._tree)
        edge.attach(self.root, node)

    # :4.2;
    #
    #   .
    #   |4.2
    #   .
    cpdef reduceB(self, Newick.TokenColon colon,
      Newick.TokenBranchLength branchLength, Newick.TokenSemicolon semicolon):
        "%accept"
        cdef Node node
        cdef Edge edge

        node = Node(self.parser._tree)
        edge = Edge(self.parser._tree)
        edge.lengthSet(float(branchLength.raw))
        edge.attach(self.root, node)

    # ;
    #
    #  .
    cpdef reduce(self, Newick.TokenSemicolon semicolon):
        "%accept"

cdef class _NewickDescendantList(Newick.DescendantList):
    "%extend DescendantList"

    cdef Node node

    cpdef reduce(self, Newick.TokenLparen lparen,
      Newick.SubtreeList SubtreeList, Newick.TokenRparen rparen):
        "%accept"
        cdef _NewickSubtree subtree
        cdef Edge edge

        self.node = Node(self.parser._tree)
        for subtree in (<_NewickSubtreeList>SubtreeList).subtrees:
            edge = Edge(self.parser._tree)
            edge.lengthSet(subtree.len)
            edge.attach(self.node, subtree.node)

cdef class _NewickSubtreeList(Newick.SubtreeList):
    "%extend SubtreeList"

    cdef list subtrees

    cpdef reduceOne(self, Newick.Subtree Subtree):
        "%accept"
        self.subtrees = [Subtree]

    cpdef reduceExtend(self, Newick.SubtreeList SubtreeList,
      Newick.TokenComma comma, Newick.Subtree Subtree):
        "%accept"
        self.subtrees = (<_NewickSubtreeList>SubtreeList).subtrees
        self.subtrees.insert(0, Subtree)

cdef class _NewickSubtree(Newick.Subtree):
    "%extend Subtree"

    cdef Node node
    cdef float len

    cpdef reduceDIB(self, Newick.DescendantList DescendantList,
      Newick.Label Label, Newick.TokenColon colon,
      Newick.TokenBranchLength branchLength):
        "%accept"
        self.node = (<_NewickDescendantList>DescendantList).node
        self.parser._labelNode(self.node, <_NewickLabel>Label)
        self.len = float(branchLength.raw)

    cpdef reduceDI(self, Newick.DescendantList DescendantList,
      Newick.Label Label):
        "%accept"
        self.node = (<_NewickDescendantList>DescendantList).node
        self.parser._labelNode(self.node, <_NewickLabel>Label)
        self.len = 0.0

    cpdef reduceDB(self, Newick.DescendantList DescendantList,
      Newick.TokenBranchLength branchLength):
        "%accept"
        self.node = (<_NewickDescendantList>DescendantList).node
        self.len = float(branchLength.raw)

    cpdef reduceLB(self, Newick.Label Label, Newick.TokenColon colon,
      Newick.TokenBranchLength branchLength):
        "%accept"
        self.node = Node(self.parser._tree)
        self.parser._labelNode(self.node, <_NewickLabel>Label)
        self.len = float(branchLength.raw)

    cpdef reduceL(self, Newick.Label Label):
        "%accept"
        self.node = Node(self.parser._tree)
        self.parser._labelNode(self.node, <_NewickLabel>Label)
        self.len = 0.0

cdef class _NewickLabel(Newick.Label):
    "%extend Label"

    cdef str raw

    cpdef reduceU(self, Newick.TokenUnquotedLabel unquotedLabel):
        "%accept"
        self.raw = unquotedLabel.raw

    cpdef reduceQ(self, Newick.TokenQuotedLabel quotedLabel):
        "%accept"
        self.raw = quotedLabel.raw

    cpdef reduceB(self, Newick.TokenBranchLength branchLength):
        "%accept"
        self.raw = branchLength.raw

    cpdef reduceE(self):
        "%accept"
        self.raw = None

cdef Parsing.Spec _NewickSpec

cdef class _NewickParser(Newick.Parser):
    cdef readonly Tree _tree
    cdef TaxonMap _taxonMap
    cdef bint _newickAutoMap

    def __init__(self, Tree tree, TaxonMap taxonMap, bint newickAutoMap=False):
        global _NewickSpec

        if _NewickSpec is None:
            _NewickSpec = self._initSpec()
        Newick.Parser.__init__(self, _NewickSpec)

        self._tree = tree
        self._taxonMap = taxonMap
        self._newickAutoMap = newickAutoMap

    cdef Parsing.Spec _initSpec(self):
        return Parsing.Spec([sys.modules[__name__], Crux.Newick],
          startSym=Newick.Tree, pickleFile="%s/share/Crux-%s/Tree.pickle" %
          (Crux.Config.prefix, Crux.Config.version),
          verbose=(False if (not __debug__ or Crux.opts.quiet) else True),
          skinny=(False if __debug__ else True),
          logFile="%s/share/Crux-%s/Tree.log" %
          (Crux.Config.prefix, Crux.Config.version))

    cpdef _labelNode(self, Node node, _NewickLabel label):
        cdef int ind

        if label.raw is None:
            return

        if self._newickAutoMap:
            ind = self._taxonMap.ntaxaGet()
            self._taxonMap.map(label.raw, ind)
        else:
            ind = self._taxonMap.indGet(label.raw)
            if ind == -1:
                raise Malformed("No TaxonMap entry for %r" % label.raw)
        node.taxonNumSet(ind)

    # Overridden method.
    cpdef parse(self, str input, int begPos=0, int line=1, int col=0,
      bint verbose=False):
        cdef ret = Newick.Parser.parse(self, input, begPos, line, col, verbose)

        return ret

#
# End Newick tree construction support.
#===============================================================================

# Default branch length callback function for random tree construction.
def _defaultRandomBranchCallback():
    return 1.0

# XXX Add some sort of taxon re-numbering method?
cdef class Tree:
    def __init__(self, with_=None, TaxonMap taxonMap=None,
      bint unrooted=True, bint newickAutoMap=False, randomBranchCallback=None):
        self._sn = 0
        self._cacheSn = -1
        self.unrooted = unrooted
        if type(with_) == int:
            self._randomNew(with_, taxonMap, randomBranchCallback)
        elif type(with_) == str:
            if taxonMap == None:
                taxonMap = TaxonMap()
            self._taxonMap = taxonMap
            self._newickNew(with_, unrooted, newickAutoMap)
        else:
            assert with_ is None
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

    cdef void _newickNew(self, str input, bint unrooted, bint newickAutoMap) \
      except *:
        cdef _NewickParser parser
        cdef _NewickTree tree

        parser = _NewickParser(self, self._taxonMap, newickAutoMap)
        parser.parse(input)
        tree = parser.start[0]
        self.baseSet(tree.root)
        if unrooted:
            self.unrooted = False
            self.deroot()

    def _dup(self, newTree, node, prevRing):
        newNode = Node(newTree)
        taxonNum = node._taxonNum
        if taxonNum != -1:
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

    # XXX Cache node degrees?
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

    cpdef deroot(self):
        cdef Node node
        cdef Edge edge
        cdef Ring ring
        cdef float removedLength

        if self.unrooted:
            return
        self.unrooted = True

        ring = self._base._ring
        if ring is not None:
            ring = ring.other()
            node = ring._node
            edge = ring._edge
            edge.detach()
            if node.degree() == 2:
                ring = node._ring
                self.baseSet(ring.other()._node)
                edge = ring._edge
                removedLength = edge._length
                edge.detach()
                ring = node._ring
                edge = ring._edge
                node = ring.other()._node
                edge.detach()
                # XXX Change ring header/order in order to disturb canonical
                # order?
                edge.attach(self._base, node)
                edge.lengthSet(edge._length + removedLength)
            else:
                self.baseSet(node)
        else:
            # Only the root node exists, so discard the whole tree.
            self.baseSet(None)

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

    cpdef render(self, bint labels=False, bint lengths=False,
      lengthFormat="%.5e", outFile=None):
        cdef ret
        cdef Node n, neighbor
        cdef Ring ring
        cdef int degree

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
            if self.unrooted:
#                callback("[&u] ")
                degree = n.degree()
                if degree == 0:
                    # There is only one node in the tree.
                    n.rrender(None, self._taxonMap, labels, lengths,
                      lengthFormat, callback)
                elif degree == 1:
                    # Leaf node.  If this node's neighbor is an internal node,
                    # start rendering with it, in order to unroot the tree.
                    ring = n.ring()
                    assert ring != None
                    neighbor = ring.other()._node
                    if neighbor.degree() > 1:
                        # Start with the internal node.
                        neighbor.rrender(None, self._taxonMap, labels,
                          lengths, lengthFormat, callback, noLength=True)
                    else:
                        # This tree only has two taxa; start with the tree
                        # base.
                        callback("(")
                        n.rrender(neighbor, self._taxonMap, labels, lengths,
                          lengthFormat, callback)
                        callback(",")
                        neighbor.rrender(n, self._taxonMap, labels, lengths,
                          lengthFormat, callback, zeroLength=True)
                        callback(")")
                else:
                    # Internal node.
                    n.rrender(None, self._taxonMap, labels, lengths,
                      lengthFormat, callback, noLength=True)
            else: # Rooted tree.
#                callback("[&r] ")
                if n._taxonNum != -1:
                    raise Malformed("Root is labeled")

                degree = n.degree()
                if degree > 1:
                    raise Malformed("Root is an internal node")

                if degree == 1:
                    ring = n.ring()
                    assert ring != None
                    neighbor = ring.other()._node
                    neighbor.rrender(n, self._taxonMap, labels, lengths,
                      lengthFormat, callback)
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

    cpdef rrender(self, Node prev, TaxonMap taxonMap, bint labels, bint lengths,
      lengthFormat, callback, bint zeroLength=False, bint noLength=False):
        cdef bint did_paren = False
        cdef object taxonLabel, m
        cdef Ring ring, r
        cdef Node neighbor
        cdef int degree

        # Iterate through neighbors.
        ring = self._ring
        if ring != None:
            for r in ring:
                # Get the node on the other end of the edge.  If it isn't prev,
                # recurse.
                neighbor = r.other()._node
                if neighbor != prev:
                    if did_paren:
                        callback(",")
                    elif not did_paren:
                        callback("(")
                        did_paren = True

                    neighbor.rrender(self, taxonMap, labels, lengths,
                      lengthFormat, callback)

            if did_paren:
                callback(")")

        # Render label.
        if self._taxonNum != -1:
            if labels:
                # Protect special characters, if necessary.
                taxonLabel = taxonMap.labelGet(self._taxonNum)
                m = re.compile(r"[_()[\]':;,]").search(taxonLabel)
                if m:
                    taxonLabel = re.compile("'").sub("''", taxonLabel)
                    callback("'%s'" % taxonLabel)
                else:
                    if taxonLabel.find(" ") != -1:
                        taxonLabel = re.compile(" ").sub("_", taxonLabel)
                    callback("%s" % taxonLabel)
            else:
                callback("%d" % self._taxonNum)

        # Render branch length.
        degree = self.degree()
        if lengths and degree > 0:
            ring = self._ring
            assert ring != None
            if zeroLength:
                # This tree only has two taxa; take care not to double the
                # branch length.
                callback((":" + lengthFormat) % 0.0)
            elif not noLength:
                callback((":" + lengthFormat) % (ring._edge._length))

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
