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
    m = crux.taxon_map.taxon_map(map)
    print
    sys.stdout.write("map: ")
    print m.taxa_get()
    print "ntaxa: %d" % m.ntaxa_get()

    i = 0
    for label in m.taxa_get():
        print "label '%s' '%s' at index %d %d" \
              % (label, m.label_get(i), m.ind_get(label), i)

        i += 1

m = crux.taxon_map.taxon_map()
print
sys.stdout.write("map: ")
print m.taxa_get()
print "ntaxa: %d" % m.ntaxa_get()

m.map('B', 1)
m.map('A', 0)
sys.stdout.write("map: ")
print m.taxa_get()
print "ntaxa: %d" % m.ntaxa_get()

m.map('D', 3)
sys.stdout.write("map: ")
try:
    print m.taxa_get()
except ValueError:
    print "ValueError in m.taxa_get()"
print "ntaxa: %d" % m.ntaxa_get()

m.map('C', 2)
sys.stdout.write("map: ")
print m.taxa_get()
print "ntaxa: %d" % m.ntaxa_get()

print "Test end"
