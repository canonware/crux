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
import newick
import taxon_map

import random

class _newick_parser(newick.newick):
    def __init__(self, tree, map):
        self._tree = tree
        self._map = map
        self._taxon_stack = []
        pass

    # Overridden method.
    def parse(self, input, tree):
        if not newick.newick.parse(self, input):
            if len(self._taxon_stack) > 0:
                tree.base_set(self._taxon_stack[0])
            retval = False
        else:
            retval = True

        return retval

    # Overridden method.
    def open_paren_accept(self):
        self._taxon_stack.insert(0, None)

    # Overridden method.
    def close_paren_accept(self):
        # Create an internal node, and join the top nodes to it.
        cnt = 0
        for elm in self._taxon_stack:
            if elm == None:
                break
            elif type(elm) == node.node:
                cnt += 1
        if cnt < 2:
            # Not enough neighbors on stack to join nodes together.  Remove open
            # paren (represented as None) from stack.
            self._taxon_stack.remove(None)
        else:
            # Create new node.
            nnode = node.node(self._tree)

            # Iteratively connect neighboring nodes to nnode.
            i = 0
            while i < cnt:
                if type(self._taxon_stack[0]) == float:
                    length = self._taxon_stack.pop(0)
                else:
                    length = 0.0

                n = self._taxon_stack.pop(0)
                e = edge.edge(self._tree)
                e.attach(nnode, n)
                e.length_set(length)
                i += 1

            # Pop paren (None).
            self._taxon_stack.pop(0)

            # Push nnode onto the stack.
            self._taxon_stack.insert(0, nnode)

    # Helper method, called by {root,leaf}_label_accept().
    def _label_accept(self):
        if self._map.ind_get(self.token()) != None:
            # Taxon mapping defined.
            val = self._map.ind_get(self.token())
        else:
            # No taxon mapping defined; try to convert the label to an
            # integer.
            try:
                val = int(self.token())
                self._map.map(self.token(), val)
            except ValueError:
                # Failed conversion.  Create a new mapping.
                self._map.map(self.token(), self._map.ntaxa_get())

        # Create a new node and push it onto the stack.
        nnode = node.node(self._tree)
        nnode.taxon_num_set(val)
        self._taxon_stack.insert(0, nnode)

    # Overridden method.
    def root_label_accept(self):
        # A tree with only one node can be written as a single root label.  This
        # is the only case that we allow root labels here.
        if len(self.token()) > 0:
            if len(self._taxon_stack) > 0:
                self.error_raise("Trailing root label not supported")
            else:
                self._label_accept()

    # Overridden method.
    def leaf_label_accept(self):
        self._label_accept()

    # Overridden method.
    def length_accept(self):
        self._taxon_stack.insert(0, float(self.token()))

    # Overridden method.
    def semicolon_accept(self):
        # If there is an internal node with only two neighbors on the stack,
        # splice it out of the tree.
        if len(self._taxon_stack) != 0:
            if type(self._taxon_stack[0]) == float:
                length = self._taxon_stack.pop(0)
            else:
                length = None

            n = self._taxon_stack[0]
            if n.degree() == 2:
                self._taxon_stack.pop(0)
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
                    e.length_set(length)
                else:
                    e.length_set(edge_a.length_get() + edge_b.length_get())
                    
                # Push a node back onto the stack.
                self._taxon_stack.insert(0, node_a)

class tree(_tree.Tree):
    def __init__(self, with=None, map=taxon_map.taxon_map()):
        self._map = map

        if type(with) == int:
            self._random_new(with)
        elif type(with) == str or type(with) == file:
            self._newick_new(with)
        elif type(with) == list:
            self._nj_new(with)

    def _random_new(self, ntaxa):
        # Create a stack of leaf nodes.
        subtrees = []
        i = 0
        while i < ntaxa:
            nnode = node.node(self)
            nnode.taxon_num_set(i)
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

        self.base_set(subtree_a)

    def _newick_new(self, input):
        parser = _newick_parser(self, self._map)
        return parser.parse(input, self)

    def _nj_new(self, input):
        self._nj(input)

    def prints(self, labels=False, lengths=False):
        n = self.base_get()
        if n != None:
            retval = n.rprints(None, self._map, labels, lengths)

            if n.taxon_num_get() == None:
                # Internal node.
                retval = "%s;" % retval
            else:
                # Leaf node.
                retval = "(%s);" % retval
        else:
            retval = ";"

        return retval
#EOF
