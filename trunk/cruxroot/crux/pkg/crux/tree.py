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

import _tree

import node
import edge
import NewickParser
import TaxonMap

import random

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
            elif type(elm) == node.node:
                cnt += 1
        if cnt < 2:
            # Not enough neighbors on stack to join nodes together.  Remove open
            # paren (represented as None) from stack.
            self._taxonStack.remove(None)
        else:
            # Create new node.
            nnode = node.node(self._tree)

            # Iteratively connect neighboring nodes to nnode.
            i = 0
            while i < cnt:
                if type(self._taxonStack[0]) == float:
                    length = self._taxonStack.pop(0)
                else:
                    length = 0.0

                n = self._taxonStack.pop(0)
                e = edge.edge(self._tree)
                e.attach(nnode, n)
                e.lengthSet(length)
                i += 1

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
        nnode = node.node(self._tree)
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
                # Get edges (and end that n is attached to).
                (edge_a, end_a) = n.edge()
                (edge_b, end_b) = edge_a.next(end_a)
                # Get neighboring nodes.
                node_a = edge_a.node(end_a ^ 1)
                node_b = edge_b.node(end_b ^ 1)
                # Detach edges.
                edge_a.detach()
                edge_b.detach()
                # Attach neighbors.
                e = edge.edge(self._tree)
                e.attach(node_a, node_b)
                if length != None:
                    e.lengthSet(length)
                else:
                    e.lengthSet(edge_a.lengthGet() + edge_b.lengthGet())
                    
                # Push a node back onto the stack.
                self._taxonStack.insert(0, node_a)

class tree(_tree.Tree):
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
        i = 0
        while i < ntaxa:
            nnode = node.node(self)
            nnode.taxonNumSet(i)
            subtrees.append(nnode)
            self._map.map(str(i), i)

            i += 1

        # Iteratively randomly remove two items from the stack, join them, and
        # push the result back onto the stack.  Stop when there are two subtrees
        # left.
        while len(subtrees) > 2:
            subtree_a = subtrees.pop(random.randint(0, len(subtrees) - 1))
            subtree_b = subtrees.pop(random.randint(0, len(subtrees) - 1))
            nnode = node.node(self)
            edge.edge(self).attach(nnode, subtree_a)
            edge.edge(self).attach(nnode, subtree_b)
            subtrees.append(nnode)

        # Attach the last two subtrees directly, in order to finish constructing
        # an unrooted tree.
        subtree_a = subtrees.pop(0)
        subtree_b = subtrees.pop(0)
        edge.edge(self).attach(subtree_a, subtree_b)

        self.baseSet(subtree_a)

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
