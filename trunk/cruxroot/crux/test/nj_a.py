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

print "Test begin"

matrix = crux.DistMatrix.DistMatrix(open('test/treezilla.dist'))
t = crux.Tree.Tree(matrix)

t.canonize()
print t.prints()

print "Test end"
