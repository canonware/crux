################################################################################
#
# <Copyright = jasone>
# <License>
#
################################################################################
#
# Version: Crux <Version = crux>
#
################################################################################

from C_Tree import *

import Node
import Edge
import Ring
import NewickParser
import TaxonMap
import crux
import __builtin__

import random

class _NewickParser(NewickParser.NewickParser):
    def __init__(self, tree, map, autoMap=False):
        self._tree = tree
        self._map = map
        self._taxonStack = []
        self._autoMap = autoMap

    # Overridden method.
    def parse(self, input):
        if not NewickParser.NewickParser.parse(self, input):
            if len(self._taxonStack) > 0:
                self._tree.baseSet(self._taxonStack[0])
            retval = False
        else:
            retval = True

        return retval

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
            elif type(elm) == Node.Node:
                cnt += 1
        if cnt < 2:
            # Not enough neighbors on stack to join nodes together.  Remove open
            # paren (represented as None) from stack.
            self._taxonStack.remove(None)
        else:
            # Create new node.
            nnode = Node.Node(self._tree)

            # Iteratively connect neighboring nodes to nnode.
            for i in forints(cnt):
                if type(self._taxonStack[0]) == float:
                    length = self._taxonStack.pop(0)
                else:
                    length = 0.0

                n = self._taxonStack.pop(0)
                e = Edge.Edge(self._tree)
                e.attach(nnode, n)
                e.lengthSet(length)

            # Pop paren (None).
            self._taxonStack.pop(0)

            # Push nnode onto the stack.
            self._taxonStack.insert(0, nnode)

    # Helper method, called by {root,leaf}LabelAccept().
    def _labelAccept(self):
        if self._map.indGet(self.token()) != None:
            # Taxon mapping defined.
            val = self._map.indGet(self.token())
        else:
            # No taxon mapping defined; try to convert the label to an
            # integer.
            try:
                val = int(self.token())
                self._map.map(self.token(), val)
            except __builtin__.ValueError:
                if self._autoMap:
                    # Create a new mapping.
                    val = self._map.ntaxaGet()
                    self._map.map(self.token(), val)
                else:
                    # Failed conversion.
                    raise crux.Tree.ValueError, "No mapping for '%s'" \
                          % self.token()

        # Create a new node and push it onto the stack.
        nnode = Node.Node(self._tree)
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
                e = Edge.Edge(self._tree)
                e.attach(nodeA, nodeB)
                if length != None:
                    e.lengthSet(length)
                else:
                    e.lengthSet(ringA.edge().lengthGet()
                                + ringB.edge().lengthGet())

                # Push a node back onto the stack.
                self._taxonStack.insert(0, nodeA)

class Tree(C_Tree):
    def __init__(self, with=None, map=None, autoMap=False):
        if type(with) == int:
            self._randomNew(with, map)
        elif type(with) == str or type(with) == file:
            if map == None:
                map = TaxonMap.TaxonMap()
            self._map = map
            self._newickNew(with, autoMap)
        else:
            if map == None:
                map = TaxonMap.TaxonMap()
            self._map = map

    def _randomNew(self, ntaxa, map):
        if map == None:
            self._map = TaxonMap.TaxonMap()
        else:
            self._map = map

        # Create a stack of leaf nodes.
        subtrees = []
        for i in forints(ntaxa):
            nnode = Node.Node(self)
            nnode.taxonNumSet(i)
            subtrees.append(nnode)
            if map == None:
                self._map.map("T%d" % i, i)

        # Iteratively randomly remove two items from the stack, join them, and
        # push the result back onto the stack.  Stop when there are two subtrees
        # left.
        while len(subtrees) > 2:
            subtreeA = subtrees.pop(random.randint(0, len(subtrees) - 1))
            subtreeB = subtrees.pop(random.randint(0, len(subtrees) - 1))
            nnode = Node.Node(self)
            edgeA = Edge.Edge(self)
            edgeA.lengthSet(random.expovariate(1))
            edgeA.attach(nnode, subtreeA)
            edgeB = Edge.Edge(self)
            edgeB.lengthSet(random.expovariate(1))
            edgeB.attach(nnode, subtreeB)
            subtrees.append(nnode)

        subtreeA = subtrees.pop(0)
        if (len(subtrees) > 0):
            # Attach the last two subtrees directly, in order to finish
            # constructing an unrooted tree.
            subtreeB = subtrees.pop(0)
            edge = Edge.Edge(self)
            edge.lengthSet(random.expovariate(1))
            edge.attach(subtreeA, subtreeB)

        self.baseSet(subtreeA)

    def _newickNew(self, input, autoMap):
        parser = _NewickParser(self, self._map, autoMap)
        return parser.parse(input)

    def taxonMapGet(self):
        return self._map

    def rf(self, other):
        if type(other) == Tree:
            return self._rfPair(other)
        else:
            return self._rfSequence(other)

    # Render the tree in Newick format to a string or a file.
    def render(self, labels=False, lengths=False, lengthFormat="%.5e",
               outFile=None):
        # Set up for sending output to either a string or a file.
        if outFile == None:
            callback = self._stringRenderCallback
            self._string = ""
        else:
            callback = self._fileRenderCallback
            self._outFile = outFile

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
                        # Start with the internal node.
                        neighbor.rrender(None, self._map, labels, lengths,
                                         lengthFormat, callback)
                        callback(";")
                    else:
                        # This tree only has two taxa; start with the tree base.
                        callback("(")
                        n.rrender(None, self._map, labels, lengths,
                                  lengthFormat, callback, twoTaxa=True)
                        callback(");")
                else:
                    # There is only one node in the tree.
                    n.rrender(None, self._map, labels, lengths, lengthFormat,
                              callback)
                    callback(";")
            else:
                # Internal node.
                n.rrender(None, self._map, labels, lengths, lengthFormat,
                          callback)
                callback(";")
        else:
            callback(";")

        # Clean up and set retval according to where the output was sent.
        if outFile == None:
            retval = self._string
            self._string = None
        else:
            callback("\n")
            retval = None
            self._outFile = None
        return retval

    # Callback method that is used by the render method for recursive rendering
    # of the tree in Newick format.
    def _stringRenderCallback(self, string):
        # Append string to previous strings that were passed to this callback.
        self._string = "%s%s" % (self._string, string)

    # Callback method that is used by the render method for recursive rendering
    # of the tree in Newick format.
    def _fileRenderCallback(self, string):
        # Print string to self._outFile.
        self._outFile.write(string)
#EOF
