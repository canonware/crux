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

import _node

import re

class node(_node.Node):
    def __init__(self, tree):
        pass

    def rprints(self, prev, map, labels, lengths):
        retval = ""
        did_something = False
        did_paren = False

        if self.taxon_num_get() != None:
            # Leaf node.
            if labels:
                # Protect special characters, if necessary.
                taxon_label = map.label_get(self.taxon_num_get())
                m = re.compile("[^ ()[\]':;,]*[ ()[\]':;,]").match(taxon_label)
                if m:
                    taxon_label = re.compile("'").sub("''", taxon_label)
                    retval = "%s'%s'" % (retval, taxon_label)
                else:
                    retval = "%s%s" % (retval, taxon_label)
            else:
                retval = "%s%d" % (retval, self.taxon_num_get())

            if lengths:
                (edge, end) = self.edge()
                if edge != None:
                    retval = "%s:%f" % (retval, edge.length_get())
            did_something = True

        # Iterate through neighbors.
        (edge, end) = self.edge()
        if edge != None:
            i = 0
            while i < self.degree():
                # Get the node on the other end of the edge.  If it isn't prev,
                # recurse.
                neighbor = edge.node(end ^ 1)
                if neighbor != prev:
                    if did_something:
                        retval = "%s," % retval
                    elif not did_paren:
                        retval = "%s(" % retval
                        did_paren = True
                        did_something = True

                    retval = "%s%s" % \
                             (retval,
                              neighbor.rprints(self, map, labels, lengths))
                    if lengths:
                        if neighbor.taxon_num_get() == None:
                            retval = "%s:%f" % (retval, edge.length_get())

                (edge, end) = edge.next(end)
                i += 1

            if did_paren:
                retval = "%s)" % retval

        return retval
#EOF
