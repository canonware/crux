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

taxonMap = crux.TaxonMap.TaxonMap(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'])

a = crux.Tree.Tree('(A,(B,C),D);', taxonMap=taxonMap)
print "%.2f" % a.rf(a)

b = crux.Tree.Tree('(A,(B,C),D);', taxonMap=taxonMap)
print "%.2f" % a.rf(b)

c = crux.Tree.Tree('(A,B,(C,D));', taxonMap=taxonMap)
print "%.2f" % a.rf(c)

d = crux.Tree.Tree('(A,B,C,D);', taxonMap=taxonMap)
print "%.2f" % a.rf(d)

print "("
r = a.rf((a, b, c, d))
for distance in r:
    print "  %.2f" % distance
print ")"

a = crux.Tree.Tree(';', taxonMap=taxonMap)
b = crux.Tree.Tree(';', taxonMap=taxonMap)
print "%.2f" % a.rf(b)

a = crux.Tree.Tree('A;', taxonMap=taxonMap)
b = crux.Tree.Tree('A;', taxonMap=taxonMap)
print "%.2f" % a.rf(b)

a = crux.Tree.Tree('(A,B);', taxonMap=taxonMap)
b = crux.Tree.Tree('(A,B);', taxonMap=taxonMap)
print "%.2f" % a.rf(b)

try:
    a = crux.Tree.Tree('A;', taxonMap=taxonMap)
    b = crux.Tree.Tree('(A,B);', taxonMap=taxonMap)
    print "%.2f" % a.rf(b)
except:
    import sys

    error = sys.exc_info()
    print "Exception %s: %s" % (error[0], error[1])

print "Test end"
