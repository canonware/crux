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

matrix = crux.CTMatrix.CTMatrix()
f = open(opts.scriptargs[1] + "/test/treezilla.fasta")
matrix.fastaFileParse(f, 'DNA')

print matrix.fastaPrints()

print "Test end"
