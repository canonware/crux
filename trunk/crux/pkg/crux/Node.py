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

    def rrender(self, prev, map, labels, lengths, lengthFormat, outFile=None,
                twoTaxa=False):
        if outFile == None:
            retval = self._stringRrender(prev, map, labels, lengths,
                                         lengthFormat, twoTaxa)
        else:
            self._fileRrender(prev, map, labels, lengths, lengthFormat, outFile,
                              twoTaxa)
            retval = None

        return retval

    def _stringRrender(self, prev, map, labels, lengths, lengthFormat, twoTaxa):
        retval = ""
        did_something = False
        did_paren = False

        if self.taxonNumGet() != None:
            # Leaf node.
            if labels:
                # Protect special characters, if necessary.
                taxonLabel = map.labelGet(self.taxonNumGet())
                m = re.compile("[^ ()[\]':;,]*[ ()[\]':;,]").match(taxonLabel)
                if m:
                    taxonLabel = re.compile("'").sub("''", taxonLabel)
                    retval = "%s'%s'" % (retval, taxonLabel)
                else:
                    retval = "%s%s" % (retval, taxonLabel)
            else:
                retval = "%s%d" % (retval, self.taxonNumGet())

            if lengths:
                ring = self.ring()
                if ring != None:
                    if twoTaxa:
                        # This tree only has two taxa; take care not to double
                        # the branch length.
                        retval = ("%s:" + lengthFormat) \
                                 % (retval, ring.edge().lengthGet() / 2)
                    else:
                        retval = ("%s:" + lengthFormat) \
                                 % (retval, ring.edge().lengthGet())
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
                              neighbor._stringRrender(self, map, labels,
                                                      lengths, lengthFormat,
                                                      twoTaxa))

                    if lengths:
                        if neighbor.taxonNumGet() == None:
                            retval = ("%s:" + lengthFormat) \
                                     % (retval, ring.edge().lengthGet())

                ring = ring.next()

            if did_paren:
                retval = "%s)" % retval

        return retval

    def _fileRrender(self, prev, map, labels, lengths, lengthFormat, outFile,
                     twoTaxa):
        did_something = False
        did_paren = False

        if self.taxonNumGet() != None:
            # Leaf node.
            if labels:
                # Protect special characters, if necessary.
                taxonLabel = map.labelGet(self.taxonNumGet())
                m = re.compile("[^ ()[\]':;,]*[ ()[\]':;,]").match(taxonLabel)
                if m:
                    taxonLabel = re.compile("'").sub("''", taxonLabel)
                    outFile.write("'%s'" % taxonLabel)
                else:
                    outFile.write("%s" % taxonLabel)
            else:
                outFile.write("%d" % self.taxonNumGet())

            if lengths:
                ring = self.ring()
                if ring != None:
                    if twoTaxa:
                        # This tree only has two taxa; take care not to double
                        # the branch length.
                        outFile.write((":" + lengthFormat) \
                                      % (ring.edge().lengthGet() / 2))
                    else:
                        outFile.write((":" + lengthFormat) \
                                      % (ring.edge().lengthGet()))
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
                        outFile.write(",")
                    elif not did_paren:
                        outFile.write("(")
                        did_paren = True
                        did_something = True

                    neighbor._fileRrender(self, map, labels, lengths,
                                          lengthFormat, outFile, twoTaxa)

                    if lengths:
                        if neighbor.taxonNumGet() == None:
                            outFile.write((":" + lengthFormat) \
                                          % ring.edge().lengthGet())

                ring = ring.next()

            if did_paren:
                outFile.write(")")
#EOF
