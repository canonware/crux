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

matrix = crux.DistMatrix.DistMatrix(open(opts.scriptargs[1]
                                         + '/test/treezilla.dist'))
t = crux.Tree.Tree(matrix)

t.canonize()
print t.prints()

print "Test end"
