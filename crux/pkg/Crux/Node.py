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

from C_Node import *

import re

class Node(C_Node):
    def __init__(self, tree):
        pass

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
#EOF
