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

testMatrix = """4
A    8  7 12
B       9 14
C         11
D
"""

matrix = crux.DistMatrix.DistMatrix(testMatrix)
t = crux.Tree.Tree(matrix, useNj=True)
t.canonize()
t.render(labels=True, lengths=True, lengthFormat="%.6f", outFile=sys.stdout)

print "Test end"
