"""
    Toolkit for molecular phylogenetic inference.

    Crux provides various components that are critical to many phylogenetic
    analyses, such as character alignments and trees.  These components may be
    used directly to implement fundamentally new analyis algorithms, but their
    main purpose is to support Crux's higher level functionality such as
    Bayesian estimation of model parameters.

    There are three main use cases for Crux, which from most basic to most
    advanced are:

    1) Canned data analysis using one of the Python scripts that are provided
       with Crux:
         a) redpoint : Bayesian Markov chain Monte Carlo estimation of model
                       posteriors.
         b) MrRogers : Relaxed neighbor joining (RNJ) and neighbor joining (NJ)
                       tree construction, based on pairwise distances among
                       taxa.
    2) Non-canned analysis using the Python scripting interface provided by the
       crux front end script.  The crux script can be used interactively or in
       batch mode.
    3) Development of new methods, using Crux as the foundation.  It is
       possible to use Crux packages and modules via the Python and/or Cython
       programming languages.  Python is adequate if performance is not a
       concern.  Cython makes it possible to write much faster code, and it is
       even possible to cleanly integrate with C code as necessary.

    The Crux package includes multiple sub-modules and sub-packages, which are
    separately documented:

    * Crux.Character  : Phylogenetic character classes for DNA and protein.
    * Crux.Config     : Crux configuration information.
    * Crux.Copying    : Crux licensing information.
    * Crux.CTMatrix   : Character-by-taxon matrices and character alignments.
    * Crux.DistMatrix : Distance matrices.
    * Crux.Fasta      : FASTA character data format parser.
    * Crux.Mc3        : Metropolis-coupled Markov chain Monte Carlo sampler.
    * Crux.Newick     : Newick tree format parser.
    * Crux.Taxa       : Taxa and taxa map classes.
    * Crux.Tree       : Tree, node, edge, and ring classes.
"""

from Cx cimport *

# Initialize libCx.
CxInit()

# Import pure Python modules.
import Config
import Copying

# Import Cython modules.
cimport Crux.Character as Character
cimport Crux.CTMatrix as CTMatrix
cimport Crux.DistMatrix as DistMatrix
cimport Crux.Fasta as Fasta
cimport Crux.Mc3 as Mc3
cimport Crux.Newick as Newick
cimport Crux.Taxa as Taxa
cimport Crux.Tree as Tree

IF @enable_mpi@:
    from mpi4py import MPI
    cimport mpi4py.mpi_c as mpi
    cdef int size

    mpi.MPI_Comm_size(mpi.MPI_COMM_WORLD, &size)
    if size > 1:
        # Make sure that all nodes are seeded the same, even if the seed isn't
        # explicitly set by an application that directly imports the Crux
        # module.
        seed(0)

cpdef threaded():
    """
        Enable thread parallelism.  Once enabled, parallelism cannot be
        disabled without restarting Crux.  (The crux front end script calls
        this function by default.)

        When thread parallelism is enabled, Crux logically divides character
        alignments into slices and uses a pool of worker threads (one thread
        per CPU) to concurrently compute site log-likelihoods for the slices.
    """
    import Config # Work around an apparent Cython bug.
    CxThreaded()
    Config.threaded = True

cpdef seed(unsigned s):
    """
        Seed all pseudo-random number generators.
    """
    IF @enable_mpi@:
        from mpi4py import MPI
        cdef int size

        mpi.MPI_Comm_size(mpi.MPI_COMM_WORLD, &size)
        if size != 1:
            mpi.MPI_Bcast(&s, 1, mpi.MPI_UNSIGNED, 0, mpi.MPI_COMM_WORLD)

    import Config # Work around an apparent Cython bug.
    Config.seed = s
    import random
    random.seed(s)
