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

import _Node

import re

class Node(_Node.Node):
    def __init__(self, tree):
        pass

    def rprints(self, prev, map, labels, lengths):
        retval = ""
        did_something = False
        did_paren = False

        if self.taxonNumGet() != None:
            # Leaf node.
            if labels:
                # Protect special characters, if necessary.
                taxon_label = map.labelGet(self.taxonNumGet())
                m = re.compile("[^ ()[\]':;,]*[ ()[\]':;,]").match(taxon_label)
                if m:
                    taxon_label = re.compile("'").sub("''", taxon_label)
                    retval = "%s'%s'" % (retval, taxon_label)
                else:
                    retval = "%s%s" % (retval, taxon_label)
            else:
                retval = "%s%d" % (retval, self.taxonNumGet())

            if lengths:
                ring = self.ring()
                if ring != None:
                    retval = "%s:%f" % (retval, ring.edge().lengthGet())
            did_something = True

        # Iterate through neighbors.
        ring = self.ring()
        if ring != None:
            for i in forints(self.degree()):
                # Get the node on the other end of the edge.  If it isn't prev,
                # recurse.
                neighbor = ring.other().node()
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
                        if neighbor.taxonNumGet() == None:
                            retval = "%s:%f" % (retval, ring.edge().lengthGet())

                ring = ring.next()

            if did_paren:
                retval = "%s)" % retval

        return retval
#EOF
