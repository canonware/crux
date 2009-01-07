"""
    Toolkit for molecular phylogenetic inference.

    XXX Description docs here.
"""

from Cx cimport *
from SFMT cimport *

# Initialize libCx.
CxInit()

# Import pure Python modules.
import Config
import Exception

# Import Cython modules.
cimport Crux.Character as Character
cimport Crux.CTMatrix as CTMatrix
cimport Crux.DistMatrix as DistMatrix
cimport Crux.Fasta as Fasta
cimport Crux.Newick as Newick
cimport Crux.Taxa as Taxa
cimport Crux.Tree as Tree

cpdef seed(unsigned s):
    import Config
    Config.seed = s
    import random
    random.seed(s)
    init_gen_rand(s)
