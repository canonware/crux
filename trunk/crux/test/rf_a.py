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

map = crux.TaxonMap.TaxonMap(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'])

a = crux.Tree.Tree('(A,(B,C),D);', map=map)
print "%.2f" % a.rf(a)

b = crux.Tree.Tree('(A,(B,C),D);', map=map)
print "%.2f" % a.rf(b)

c = crux.Tree.Tree('(A,B,(C,D));', map=map)
print "%.2f" % a.rf(c)

print "("
r = a.rf((a, b, c))
for distance in r:
    print "  %.2f" % distance
print ")"

d = crux.Tree.Tree('A;', map=map)
e = crux.Tree.Tree('(A,B);', map=map)
print "%.2f" % d.rf(e)

print "Test end"