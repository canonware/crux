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

t = crux.Tree.Tree("(A,('B B','C''C'));",
                   crux.TaxonMap.TaxonMap(['A', 'B B', "C'C"]))
t.canonize()
print t.prints(labels=True)

print "Test end"
