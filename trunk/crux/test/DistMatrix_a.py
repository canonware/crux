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
import os
import tempfile

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
""",

"""2
Taxon_A 0.0 1.0
Taxon_B 1.0 0.0
""",

"""2
Taxon_A     1.0
Taxon_B
""",

"""2
Taxon_A
Taxon_B 1.0
""",

"""3
Taxon_A 0.0 1.0 2.0
Taxon_B 1.0 0.0 3.0
Taxon_C 2.0 3.0 0.0
""",

"""3
Taxon_A     1.0 2.0
Taxon_B         3.0
Taxon_C
""",

"""3
Taxon_A
Taxon_B 1.0
Taxon_C 2.0 3.0
""",

"""3
Taxon_A 1.0e4 +1.0e4 -1.0e4
Taxon_B 1e4   +1e4   -1e4
Taxon_C 1e+4  1e-4   1E4
""",

#
# Error cases.
#

"""
Alpha      0.000 1.000 2.000 3.000 3.000
Beta       1.000 0.000 2.000 3.000 3.000
Gamma      2.000 2.000 0.000 3.000 3.000
Delta      3.000 3.000 0.000 0.000 1.000
Epsilon    3.000 3.000 3.000 1.000 0.000
""",

"""5
           0.000 1.000 2.000 3.000 3.000
Beta       1.000 0.000 2.000 3.000 3.000
Gamma      2.000 2.000 0.000 3.000 3.000
Delta      3.000 3.000 0.000 0.000 1.000
Epsilon    3.000 3.000 3.000 1.000 0.000
""",

"""5
Alpha            1.000 2.000 3.000 3.000
Beta       1.000 0.000 2.000 3.000 3.000
Gamma      2.000 2.000 0.000 3.000 3.000
Delta      3.000 3.000 0.000 0.000 1.000
Epsilon    3.000 3.000 3.000 1.000 0.000
""",

"""5
Alpha      0.000 1.000 2.000 3.000
Beta       1.000 0.000 2.000 3.000 3.000
Gamma      2.000 2.000 0.000 3.000 3.000
Delta      3.000 3.000 0.000 0.000 1.000
Epsilon    3.000 3.000 3.000 1.000 0.000
""",

"""5
Alpha      0.000 1.000 2.000 3.000 3.000
           1.000 0.000 2.000 3.000 3.000
Gamma      2.000 2.000 0.000 3.000 3.000
Delta      3.000 3.000 0.000 0.000 1.000
Epsilon    3.000 3.000 3.000 1.000 0.000
""",

"""5
Alpha      0.000 1.000 2.000 3.000 3.000
Beta             0.000 2.000 3.000 3.000
Gamma      2.000 2.000 0.000 3.000 3.000
Delta      3.000 3.000 0.000 0.000 1.000
Epsilon    3.000 3.000 3.000 1.000 0.000
""",

"""5
Alpha      0.000 1.000 2.000 3.000 3.000
Beta       1.000 0.000 2.000 3.000 
Gamma      2.000 2.000 0.000 3.000 3.000
Delta      3.000 3.000 0.000 0.000 1.000
Epsilon    3.000 3.000 3.000 1.000 0.000
""",

"""5
Alpha      0.000 1.000 2.000 3.000 3.000
Beta       1.000 0.000
2.000 3.000 
Gamma      2.000 2.000 0.000 3.000 3.000
Delta      3.000 3.000 0.000 0.000 1.000
Epsilon    3.000 3.000 3.000 1.000 0.000
""",

"""5
Alpha      0.000 1.000 2.000 3.000 3.000
Beta       1.000 0.000 2.000 3.000 3.000
Gamma      2.000 2.000 0.000 3.000 3.000
Delta      3.000 3.000 0.000 0.000 1.000
           3.000 3.000 3.000 1.000 0.000
""",

"""5
Alpha      0.000 1.000 2.000 3.000 3.000
Beta       1.000 0.000 2.000 3.000 3.000
Gamma      2.000 2.000 0.000 3.000 3.000
Delta      3.000 3.000 0.000 0.000 1.000
Epsilon    3.000 3.000 3.000 1.000 
""",

"""5
Alpha      0.000 1.000 2.000 3.000 3.000
Beta       1.000 0.000 2.000 3.000 3.000
Gamma      2.000 2.000 0.000 3.000 3.000
Delta      3.000 3.000 0.000 0.000 1.000
Gamma      3.000 3.000 3.000 1.000 0.000
""",

            ]

print "Test begin"

for matrix in matrices:
    try:
        print "==="
        print matrix
        distMatrix = crux.DistMatrix.DistMatrix(matrix)
        print distMatrix.render('full', '%.5f')
        print distMatrix.render('upper', '%.5f')
        print distMatrix.render('lower', '%.5f')
    except:
        import sys

        error = sys.exc_info()
        print "Exception %s: %s" % (error[0], error[1])

for matrix in matrices:
    try:
        print "==="
        print matrix
        f = tempfile.TemporaryFile()
        f.write(matrix)
        f.seek(0, 0)
        distMatrix = crux.DistMatrix.DistMatrix(f)
        print distMatrix.render('full', '%.5f')
    except:
        import sys

        error = sys.exc_info()
        print "Exception %s: %s" % (error[0], error[1])

try:
    print "==="
    taxonMap = crux.TaxonMap.TaxonMap(['A', 'B', 'C'])
    distMatrix = crux.DistMatrix.DistMatrix(taxonMap)
    print distMatrix.render('full', '%.5f')
except:
    import sys

    error = sys.exc_info()
    print "Exception %s: %s" % (error[0], error[1])

print "Test end"
