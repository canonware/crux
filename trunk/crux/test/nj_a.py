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
t = matrix.nj()

t.canonize()
t.render(outFile=sys.stdout)

print "Test end"
