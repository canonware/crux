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
t = crux.Tree.Tree(matrix)
t.canonize()
print t.render(labels=True, lengths=True)

print "Test end"
