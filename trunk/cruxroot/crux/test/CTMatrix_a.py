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

dnaChar = crux.CharacterType.DnaCharacterType()
map = crux.TaxonMap.TaxonMap(['Taxon_A', 'Taxon_B', 'Taxon_C'])
matrix = crux.CTMatrix.CTMatrix(map)

matrix.charsAppend([dnaChar] * 5)
#print matrix.charsGet()
print "Character 0 (0..%d) type: %r" \
      % (len(matrix.charsGet()) - 1, type(matrix.charsGet()[0]))

print "Taxon map: %r" % matrix.taxonMapGet().taxaGet()

matrix.dataSet("Taxon_B", "CCCCC")
print "Taxon_B data: %r" % matrix.dataGet("Taxon_B")

matrix.dataSet("Taxon_A", "AAAAA")

try:
    matrix.dataSet("Taxon_C", "GGGG")
except crux.CTMatrix.Exception, x:
    import sys

    print "Exception %s: %s" % (sys.exc_type, x.__str__())

try:
    matrix.dataSet("Taxon_D", "TTTTT")
except crux.CTMatrix.Exception, x:
    import sys

    print "Exception %s: %s" % (sys.exc_type, x.__str__())

print "Test end"
