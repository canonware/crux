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

ctm = crux.CTMatrix.CTMatrix()
ctm.fastaParse("""
>a
C
>b
A
>c
C
>d
A
>e
G
""")

taxonMap = crux.TaxonMap.TaxonMap(['a', 'b', 'c', 'd', 'e'])
t = crux.Tree.Tree('((a,b),(c,(d,e)));', taxonMap=taxonMap)

t.mpPrepare(ctm)

print "Score: %d" % t.mp()

print "Test end"
