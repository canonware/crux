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

t = crux.tree("(A,('B B','C''C'));",
              crux.taxon_map.taxon_map(['A', 'B B', "C'C"]))
t.canonize()
print t.prints(labels=True)

print "Test end"
