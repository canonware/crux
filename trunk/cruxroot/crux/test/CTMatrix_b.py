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

matrix = crux.CTMatrix.CTMatrix()
f = open(opts.scriptargs[1] + "/test/treezilla.fasta")
matrix.fastaFileParse(f, 'DNA')

matrix.fastaPrint(sys.stdout)

print "Test end"
