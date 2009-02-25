"""
    Toolkit for molecular phylogenetic inference.

    XXX Description docs here.

    Crux uses parallel algorithms for some computations, such as tree
    log-likelihood.  Thread parallelism is configured by default to utilize all
    CPUs, but it is possible to disable parallelism in cases where the system
    is simultaneously running other computationally intensive applications.
"""

from Cx cimport *

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
cimport Crux.Mc3 as Mc3
cimport Crux.Newick as Newick
cimport Crux.Taxa as Taxa
cimport Crux.Tree as Tree

cpdef threaded():
    """
        Enable thread parallelism.  Once enabled, parallelism cannot be
        disabled.  (The crux front end script calls this function by default.)
    """
    import Config
    CxThreaded()
    Config.threaded = True

cpdef seed(unsigned s):
    """
        Seed all pseudo-random number generators.
    """
    import Config
    Config.seed = s
    import random
    random.seed(s)
