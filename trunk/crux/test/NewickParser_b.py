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

t = crux.Tree.Tree("(A,('B B','C''C'));",
                   crux.TaxonMap.TaxonMap(['A', 'B B', "C'C"]))
t.canonize()
print t.render(labels=True)
t.render(labels=True, outFile=sys.stdout)

print "Test end"
