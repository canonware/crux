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

testMatrices = (
"""8
A     4 13 13 10 14 25 23
B       15 15 12 16 27 25
C          16  7 11 22 20
D             13 17 28 26
E                 6 17 15
F                   17 15
G                      10
H
""",

"""8
A  0  4 13 13 10 14 25 23
B  4  0 15 15 12 16 27 25
C 13 15  0 16  7 11 22 20
D 13 15 16  0 13 17 28 26
E 10 12  7 13  0  6 17 15
F 14 16 11 17  6  0 17 15
G 25 27 22 28 17 17  0 10
H 23 25 20 26 15 15 10  0
"""
)

for testMatrix in testMatrices:
    matrix = crux.DistMatrix.DistMatrix(testMatrix)
    #matrix.render(distFormat="%2.0f", outFile=sys.stdout)
    t = crux.Tree.Tree(matrix)
    map = t.taxonMapGet()
    labels = map.taxaGet()
    t.canonize()
    t.render(labels=True, lengths=True, lengthFormat="%.0f", outFile=sys.stdout)

    for i in forints(10):
        matrix = crux.DistMatrix.DistMatrix(testMatrix)
        matrix.shuffle()
        #matrix.render(distFormat="%2.0f", outFile=sys.stdout)
        t = crux.Tree.Tree(matrix)

        treeStr = t.render(labels=True, lengths=True, lengthFormat="%.0f")
        t = crux.Tree.Tree(treeStr, map=map)
        t.canonize()
        t.render(labels=True, lengths=True, lengthFormat="%.0f",
                 outFile=sys.stdout)

print "Test end"
