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
A  0  8  7 12
B  8  0  9 14
C  7  9  0 11
D 12 14 11  0
"""

matrix = crux.DistMatrix.DistMatrix(testMatrix)
t = crux.Tree.Tree(matrix)
t.canonize()
t.render(labels=True, lengths=True, lengthFormat="%.6f", outFile=sys.stdout)
print t.render(labels=True, lengths=True, lengthFormat="%.6f")

print "Test end"
