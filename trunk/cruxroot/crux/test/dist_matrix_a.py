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

def matrix_print(matrix, ntaxa):
    x = 0
    while x < ntaxa:
        y = 0
        while y < ntaxa:
            sys.stdout.write(" %1.3f" % matrix[x * ntaxa + y])
            y += 1
        sys.stdout.write("\n")
        x += 1

matrices = ["""     5
Alpha      0.000 1.000 2.000 3.000 3.000
Beta       1.000 0.000 2.000 3.000 3.000
Gamma      2.000 2.000 0.000 3.000 3.000
Delta      3.000 3.000 0.000 0.000 1.000
Epsilon    3.000 3.000 3.000 1.000 0.000
""",

"""5
Taxon_A 0.0 1.0 2.0 3.0 4.0
Taxon_B 1.0 0.0 1.5 2.5 3.5
Taxon_C 2.0 1.5 0.0 2.2 3.2
Taxon_D 3.0 2.5 2.2 0.0 3.1
Taxon_E 4.0 3.5 3.2 3.1 0.0
""",

"""5
Taxon_A     1.0 2.0 3.0 4.0
Taxon_B         1.5 2.5 3.5
Taxon_C             2.2 3.2
Taxon_D                 3.1
Taxon_E
""",

"""5
Taxon_A
Taxon_B 1.0
Taxon_C 2.0 1.5
Taxon_D 3.0 2.5 2.2
Taxon_E 4.0 3.5 3.2 3.1
"""]

print "Test begin"

for matrix in matrices:
    (map, matrix) = crux.dist_matrix().parse(matrix)

    print map.taxa_get()
    matrix_print(matrix, map.ntaxa_get())

print "Test end"
