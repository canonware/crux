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

maps = [
    [],
    ['A'],
    ['A', 'B', 'C', 'D']
    ]

for map in maps:
    m = crux.TaxonMap.TaxonMap(map)
    print
    print "map:", m.taxaGet()
    print "ntaxa: %d" % m.ntaxaGet()

    i = 0
    for label in m.taxaGet():
        print "label '%s' '%s' at index %d %d" \
              % (label, m.labelGet(i), m.indGet(label), i)

        i += 1

m = crux.TaxonMap.TaxonMap()
print
print "map:", m.taxaGet()
print "ntaxa: %d" % m.ntaxaGet()

m.map('B', 1)
m.map('A', 0)
print "map:", m.taxaGet()
print "ntaxa: %d" % m.ntaxaGet()

m.map('D', 3)
try:
    print "map:", m.taxaGet()
except:
    import sys

    error = sys.exc_info()
    print "Exception %s: %s" % (error[0], error[1])
print "ntaxa: %d" % m.ntaxaGet()

m.map('C', 2)
print "map:", m.taxaGet()
print "ntaxa: %d" % m.ntaxaGet()

print "Test end"
