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

import sys

print "Test begin"

matrix = crux.DistMatrix.DistMatrix(open(opts.scriptargs[1]
                                         + '/test/treezilla.dist'))
t = crux.Tree.Tree(matrix, useNj=True)

t.canonize()
t.render(outFile=sys.stdout)

print "Test end"
