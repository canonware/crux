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

import _Tree

import Node
import Edge
import Ring
import NewickParser
import TaxonMap

import random

class Exception(_Tree.Exception):
    pass

class ValueError(Exception, _Tree.ValueError):
    pass

class _NewickParser(NewickParser.NewickParser):
    def __init__(self, tree, map):
        self._tree = tree
        self._map = map
        self._taxonStack = []
        pass

    # Overridden method.
    def parse(self, input, tree):
        if not NewickParser.NewickParser.parse(self, input):
            if len(self._taxonStack) > 0:
                tree.baseSet(self._taxonStack[0])
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
            except ValueError:
                # Failed conversion.  Create a new mapping.
                # XXX This allows horrible failures to easily happen later on,
                # and it's difficult to figure out what went wrong.  Is this a
                # good idea?  Perhaps raising an exception is a better way to
                # go.
                print "Create new mapping for '%s' (dangerous)" % self.token()
                val = self._map.ntaxaGet()
                self._map.map(self.token(), val)

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

class Tree(_Tree.Tree):
    def __init__(self, with=None, map=TaxonMap.TaxonMap()):
        self._map = map

        if type(with) == int:
            self._randomNew(with)
        elif type(with) == str or type(with) == file:
            self._newickNew(with)
        elif type(with) == list:
            self._njNew(with)

    def _randomNew(self, ntaxa):
        # Create a stack of leaf nodes.
        subtrees = []
        for i in forints(ntaxa):
            nnode = Node.Node(self)
            nnode.taxonNumSet(i)
            subtrees.append(nnode)
            self._map.map(str(i), i)

        # Iteratively randomly remove two items from the stack, join them, and
        # push the result back onto the stack.  Stop when there are two subtrees
        # left.
        while len(subtrees) > 2:
            subtreeA = subtrees.pop(random.randint(0, len(subtrees) - 1))
            subtreeB = subtrees.pop(random.randint(0, len(subtrees) - 1))
            nnode = Node.Node(self)
            Edge.Edge(self).attach(nnode, subtreeA)
            Edge.Edge(self).attach(nnode, subtreeB)
            subtrees.append(nnode)

        # Attach the last two subtrees directly, in order to finish constructing
        # an unrooted tree.
        subtreeA = subtrees.pop(0)
        subtreeB = subtrees.pop(0)
        Edge.Edge(self).attach(subtreeA, subtreeB)

        self.baseSet(subtreeA)

    def _newickNew(self, input):
        parser = _NewickParser(self, self._map)
        return parser.parse(input, self)

    def _njNew(self, input):
        self._nj(input)

    def prints(self, labels=False, lengths=False):
        n = self.baseGet()
        if n != None:
            retval = n.rprints(None, self._map, labels, lengths)

            if n.taxonNumGet() == None:
                # Internal node.
                retval = "%s;" % retval
            else:
                # Leaf node.
                retval = "(%s);" % retval
        else:
            retval = ";"

        return retval
#EOF
