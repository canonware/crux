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

t = crux.tree.tree("(A,('B B','C''C'));",
                   crux.TaxonMap.TaxonMap(['A', 'B B', "C'C"]))
t.canonize()
print t.prints(labels=True)

print "Test end"
