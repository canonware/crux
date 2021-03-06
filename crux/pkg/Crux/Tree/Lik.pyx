"""
    Models of molecular evolution and tree likelihoods.
"""
import cPickle
import random

import Crux.Config as Config
from Crux.Character cimport Character
from Crux.Taxa cimport Taxon
from Crux.Tree cimport Tree, Node, Edge, Ring
from Crux.CTMatrix cimport Alignment

from Cx cimport CxNcpus
from libc cimport *
from libm cimport *
from CxLik cimport *
from CxMath cimport *

IF @enable_mpi@:
    from mpi4py import MPI
    cimport mpi4py.mpi_c as mpi

DEF LikDebug = False

# Conditional likelihoods are allocated such that they are aligned to cacheline
# boundaries, and striping is configured to avoid false cacheline sharing among
# worker threads.  Over-estimating cacheline size is okay (within reason), as
# long as the estimate is a multiple of the true cacheline size.
cdef enum:
    cacheLine = 64

cdef class CL:
    """
        Conditional likelihood, associated with a ring.
    """
    def __cinit__(self):
        cdef unsigned i

        for 0 <= i < 2:
            self.cLs[i].cLMat = NULL
            self.cLs[i].lnScale = NULL
            self.cLs[i].valid = False
            self.cLs[i].parent = NULL
            self.cLs[i].nSibs = 0

    def __dealloc__(self):
        cdef unsigned i

        for 0 <= i < 2:
            if self.cLs[i].cLMat != NULL:
                free(self.cLs[i].cLMat)
                self.cLs[i].cLMat = NULL
            if self.cLs[i].lnScale != NULL:
                free(self.cLs[i].lnScale)
                self.cLs[i].lnScale = NULL

    def __init__(self):
        pass

    cdef void prepare(self, unsigned polarity, unsigned nchars, unsigned dim, \
      unsigned ncomp) except *:
        assert polarity < 2

        if self.cLs[polarity].cLMat == NULL:
            IF @have_posix_memalign@:
                if posix_memalign(<void **>&self.cLs[polarity].cLMat, \
                  cacheLine, nchars * dim * ncomp * sizeof(double)):
                    raise MemoryError("Error allocating cLMat")
            ELSE:
                self.cLs[polarity].cLMat = \
                  <double *>malloc(nchars * dim * ncomp * sizeof(double))
                if self.cLs[polarity].cLMat == NULL:
                    raise MemoryError("Error allocating cLMat")

        if self.cLs[polarity].lnScale == NULL:
            IF @have_posix_memalign@:
                if posix_memalign(<void **>&self.cLs[polarity].lnScale, \
                  cacheLine, nchars * sizeof(double)):
                    raise MemoryError("Error allocating lnScale")
            ELSE:
                self.cLs[polarity].lnScale = \
                  <double *>malloc(nchars * sizeof(double))
                if self.cLs[polarity].lnScale == NULL:
                    raise MemoryError("Error allocating lnScale")

    cdef void resize(self, unsigned polarity, unsigned nchars, unsigned dim, \
      unsigned ncomp) except *:
        cdef double *cLMat

        assert polarity < 2

        # Always (re)allocate cLMat, since ncomp may have changed.
        if self.cLs[polarity].cLMat != NULL:
            free(self.cLs[polarity].cLMat)
            self.cLs[polarity].cLMat = NULL
        self.prepare(polarity, nchars, dim, ncomp)

    cdef void flush(self, unsigned polarity) except *:
        assert polarity < 2

        if self.cLs[polarity].cLMat != NULL:
            free(self.cLs[polarity].cLMat)
            self.cLs[polarity].cLMat = NULL
        if self.cLs[polarity].lnScale != NULL:
            free(self.cLs[polarity].lnScale)
            self.cLs[polarity].lnScale = NULL
        self.cLs[polarity].valid = False
        self.cLs[polarity].parent = NULL

cdef class Lik:
    """
        Lik manages log-likelihood computation using a tree, an alignment, and
        a weighted mixture of one or more models of evolution.  Initially, the
        models are equally weighted.

        Each model of evolution consists of a Q matrix and a Gamma-distributed
        mutation rates parameter (alpha).  Q is conceptually a matrix product
        of the following form:

            Q = R * Pi
                 _            _     _            _
                |  -  a  b  c  |   | pA  0  0  0  |
                |              |   |              |
                |  a  -  d  e  |   |  0 pC  0  0  |
              = |              | * |              |,
                |  b  d  -  f  |   |  0  0 pG  0  |
                |              |   |              |
                |_ c  e  f  - _|   |_ 0  0  0 pT _|

        where {a,b,c,d,e,f} are relative substitution rates, and {pA,pC,pG,pT}
        are state frequencies.  In the general time-reversible (GTR) model,
        {a,b,c,d,e,f} can vary independently, but various sub-models can be
        chosen by specifying a constrained rate class (rclass).

        During Q matrix computation, the following conceptual steps are
        performed (though R and Pi are not directly modified):
        * A normalization factor for R is incorporated in order to adjust the
          mean instantaneous mutation rate to 1.  This affects a standardized
          interpretation of branch lengths.  (Note that the 'rmult' parameter
          is taken into consideration during normalization.  The value of
          'rmult' is irrelevant except in the context of mixture models, in
          which case it can be used to scale an entire R matrix, relative to
          other models in the mixture.)
        * A normalization factor for the state frequencies is incorporated in
          order to adjust the sum to 1.
        * The elements on Q's diagonal are set such that each row sums to 0.

        By default:
        * The rclass has only one rate ([0,0,0,0,0,0] for DNA).  The single
          rate is 1.0.
        * The state frequencies are all equal (0.25 for DNA).
        * Gamma-distributed rates are disabled (alpha is infinite and ncat is
          1).  Discretization uses category means by default, but medians can
          be used instead by setting catMedian=True when adding to the mixture
          via the constructor or addModel().

        When computing substitution probabilities for a branch of a particular
        length, numerical methods are used to compute Q's eigen decomposition,
        which can then be used to compute:

                  Q*v
        P(Q,v) = e

        For further information on the methods used, see:

          Pagel, M., A. Meade (2004) A phylogenetic mixture model for detecting
          pattern-heterogeneity in gene sequence or character-state data.
          Syst. Biol. 53(4):571-581.

          Swofford, D.L., G.J. Olsen, P.J. Waddel, D.M. Hillis (1996)
          Phylogenetic inference.  pp 407-514 in Molecular Systematics, 2nd Ed.
          (D.M. Hillis, B.K.Mable, C. Moritz, eds.).  Sinauer, Sunderland, MA.

          Yang, Z. (1994) Estimating the pattern of nucleotide substitution.
          J. Mol. Evol. 39:105-111.
    """

    def __cinit__(self):
        self.lik = NULL

    def __dealloc__(self):
        cdef CxtLik *lik
        cdef CxtLikModel *modelP
        cdef int i

        if self.lik != NULL:
            lik = self.lik
            # Iterate downward to avoid gratuitous memory moves.
            for lik.modelsLen > i >= 0:
                modelP = lik.models[i]
                self._deallocModel(modelP, i)
            free(lik.models)
            free(lik.comps)
            free(lik.siteLnL)
            free(lik.stripeLnL)
            free(lik.steps)
            free(lik)
            self.lik = NULL

    def __init__(self, Tree tree=None, Alignment alignment=None, \
      unsigned nmodels=1, unsigned ncat=1, bint catMedian=False, \
      bint invar=False):
        cdef unsigned i, nchars, npad
        cdef Character char_

        if alignment is not None:
            self._init0(tree)
            # Pad the alignment if its width isn't a multiple of the stripe
            # width.
            nchars = alignment.nchars
            npad = self._computeNpad(nchars, self._computeStripeWidth(nchars))
            char_ = alignment.charType.get()
            if npad != 0:
                alignment.pad(char_.val2code(char_.any), npad)
                nchars += npad
            self._init1(tree, nchars, char_.nstates, 0)
            self._init2(alignment, char_)
            for 0 <= i < nmodels:
                self.addModel(1.0, ncat, catMedian, invar)

    cdef unsigned _computeStripeWidth(self, unsigned nchars):
        cdef unsigned stripeQuantum, stripeWidth, nstripes

        if Config.threaded:
            # Stripe width should always be a multiple of the cacheline size in
            # order to avoid false cacheline sharing between threads.
            stripeQuantum = cacheLine / (4 * sizeof(double))
            stripeWidth = stripeQuantum

            # Limit the number of stripes.
            nstripes = nchars / stripeWidth + \
              (1 if nchars % stripeWidth != 0 else 0)
            while nstripes > CxNcpus * CxmLikMqMult:
                stripeWidth += stripeQuantum
                nstripes = nchars / stripeWidth + \
                  (1 if nchars % stripeWidth != 0 else 0)
        else:
            stripeWidth = nchars

        return stripeWidth

    cdef unsigned _computeNpad(self, unsigned nchars, unsigned stripeWidth):
        cdef unsigned npad

        if nchars % stripeWidth != 0:
            npad = stripeWidth - (nchars % stripeWidth)
        else:
            npad = 0

        return npad

    cdef void _init0(self, Tree tree) except *:
        if tree.rooted:
            raise ValueError("Tree must be unrooted")

        # Make sure there's no left over cruft from some non-likelihood
        # analyis.
        tree.clearAux()

    cdef void _init1(self, Tree tree, unsigned nchars, unsigned dim, \
      unsigned polarity) except *:
        self.mate = None
        self.tree = tree

        self.lik = <CxtLik *>calloc(1, sizeof(CxtLik))
        if self.lik == NULL:
            raise MemoryError("Error allocating lik")

        self.lik.polarity = polarity
        self.lik.dim = dim
        self.lik.rlen = self.lik.dim * (self.lik.dim-1) / 2
        self.lik.cbase = 0
        self.lik.nchars = nchars
        self.lik.mschars = nchars
        self.lik.stripeWidth = self._computeStripeWidth(nchars)
        self.lik.nstripes = self.lik.mschars / self.lik.stripeWidth
        assert self.lik.nstripes * self.lik.stripeWidth == self.lik.mschars

        IF @enable_mpi@:
            self.lik.mpiComm = mpi.MPI_COMM_NULL

        self.lik.invalidate = False
        self.lik.resize = False
        self.lik.reweight = False

        self.lik.models = <CxtLikModel **>calloc(1, sizeof(CxtLikModel *))
        if self.lik.models == NULL:
            raise MemoryError("Error allocating models")
        self.lik.modelsLen = 0
        self.lik.modelsMax = 1

        self.lik.comps = <CxtLikComp *>malloc(sizeof(CxtLikComp))
        if self.lik.comps == NULL:
            raise MemoryError("Error allocating comps")
        self.lik.compsLen = 0
        self.lik.compsMax = 1

        self.lik.siteLnL = <double *>malloc(self.lik.nchars * sizeof(double))
        if self.lik.siteLnL == NULL:
            raise MemoryError("Error allocating siteLnL")

        self.lik.stripeLnL = <double *>malloc(self.lik.nstripes * \
          sizeof(double))
        if self.lik.stripeLnL == NULL:
            raise MemoryError("Error allocating stripeLnL")

    cdef void _init2(self, Alignment alignment, Character char_) except *:
        cdef unsigned stepsMax, i

        self.char_ = char_
        self.alignment = alignment

        self.lik.npad = self.alignment.npad
        self.lik.charFreqs = self.alignment.freqs

        stepsMax = (2 * self.alignment.ntaxa) - 2
        self.lik.steps = <CxtLikStep *>malloc(stepsMax * sizeof(CxtLikStep))
        if self.lik.steps == NULL:
            raise MemoryError("Error allocating steps")
        self.lik.stepsLen = 0
        self.lik.stepsMax = stepsMax

        self.rootCL = CL()
        self.lik.rootCLC = &self.rootCL.cLs[self.lik.polarity]

    def __reduce__(self):
        return (type(self), (), self.__getstate__())

    def __getstate__(self):
        cdef initArgs, mTuple
        cdef unsigned i, j, ncat
        cdef list models, rclass, rates, freqs
        cdef double weight, rmult, alpha, wVar, wInvar
        cdef bint catMedian, invar

        initArgs = (self.tree, self.lik.nchars, self.lik.dim)
        models = []
        for 0 <= i < self.lik.modelsLen:
            weight = self.getWeight(i)
            rmult = self.getRmult(i)
            rclass = self.getRclass(i)
            rates = [self.getRate(i, j) for j in xrange(self.getNrates(i))]
            freqs = [self.getFreq(i, j) for j in xrange(self.lik.dim)]
            alpha = self.getAlpha(i)
            ncat = self.getNcat(i)
            catMedian = self.getCatMedian(i)
            invar = self.getInvar(i)
            wVar = self.getWVar(i)
            wInvar = self.getWInvar(i)

            mTuple = (weight, rmult, rclass, rates, freqs, alpha, ncat, \
              catMedian, invar, wVar, wInvar)
            models.append(mTuple)

        return (initArgs, models)

    def __setstate__(self, data):
        cdef initArgs, mTuple
        cdef unsigned i, j, m, ncat
        cdef list models, rclass, rates, freqs
        cdef double weight, rmult, alpha, wVar, wInvar
        cdef bint catMedian, invar

        (initArgs, models) = data
        self._init0(initArgs[0])
        self._init1(initArgs[0], initArgs[1], initArgs[2], 0)

        for 0 <= i < len(models):
            mTuple = models[i]
            (weight, rmult, rclass, rates, freqs, alpha, ncat, catMedian, \
              invar, wVar, wInvar) = mTuple

            m = self.addModel(weight, ncat, catMedian, invar)
            self.setRmult(i, rmult)
            self.setRclass(i, rclass, rates)

            assert len(freqs) == self.lik.dim
            for 0 <= j < self.lik.dim:
                self.setFreq(i, j, freqs[j])

            if alpha != INFINITY:
                self.setAlpha(i, alpha)

            if invar:
                self.setWVar(i, wVar)
                self.setWInvar(i, wInvar)

    cdef CxtLikModel *_allocModel(self, unsigned ncat, bint invar) except *:
        cdef CxtLikModel *modelP
        cdef CxtLikComp *comps
        cdef unsigned ncomp

        modelP = <CxtLikModel *>malloc(sizeof(CxtLikModel))
        if modelP == NULL:
            raise MemoryError("Error allocating model")

        modelP.rclass = <unsigned *>malloc(self.lik.rlen * sizeof(unsigned))
        if modelP.rclass == NULL:
            raise MemoryError("Error allocating rclass")

        modelP.rTri = <double *>malloc(self.lik.rlen * sizeof(double))
        if modelP.rTri == NULL:
            raise MemoryError("Error allocating rTri")

        modelP.piDiag = <double *>malloc(self.lik.dim * sizeof(double))
        if modelP.piDiag == NULL:
            raise MemoryError("Error allocating piDiag")

        modelP.piDiagNorm = <double *>malloc(self.lik.dim * sizeof(double))
        if modelP.piDiagNorm == NULL:
            raise MemoryError("Error allocating piDiagNorm")

        modelP.qEigVecCube = <double *>malloc(self.lik.dim * self.lik.dim * \
          self.lik.dim * sizeof(double))
        if modelP.qEigVecCube == NULL:
            raise MemoryError("Error allocating qEigVecCube")

        modelP.qEigVals = <double *>malloc(self.lik.dim * sizeof(double))
        if modelP.qEigVals == NULL:
            raise MemoryError("Error allocating qEigVals")

        ncomp = ncat
        if invar:
            ncomp += 1
        if self.lik.compsLen + ncomp > self.lik.compsMax:
            comps = <CxtLikComp *>realloc(self.lik.comps, \
              (self.lik.compsLen + ncomp) * sizeof(CxtLikComp))
            if comps == NULL:
                raise MemoryError("Error reallocating comps")
            self.lik.comps = comps
            self.lik.compsMax = (self.lik.compsLen + ncomp)
        modelP.comp0 = self.lik.compsLen
        self.lik.compsLen += ncomp
        modelP.clen = ncomp

        return modelP

    cdef void _initModel(self, CxtLikModel *modelP, double weight, \
      bint catMedian, bint invar):
        cdef unsigned i

        modelP.decomp = True
        modelP.weight = weight

        # Initialize to default equal relative rates.
        modelP.rmult = 1.0
        for 0 <= i < self.lik.rlen:
            modelP.rclass[i] = 0
            modelP.rTri[i] = 1.0

        # Initialize to default equal frequencies.
        for 0 <= i < self.lik.dim:
            modelP.piDiag[i] = 1.0 / <double>self.lik.dim

        modelP.alpha = INFINITY
        modelP.catMedian = catMedian
        modelP.invar = invar

        # Initialize comps.  Only the first one is used, unless alpha is
        # changed to something other than INFINITY, or the relative weight for
        # invariable sites is set to something other than 0.0.
        self.lik.comps[modelP.comp0].model = modelP
        self.lik.comps[modelP.comp0].cweight = 1.0
        self.lik.comps[modelP.comp0].cmult = 1.0
        for 1 <= i < modelP.clen:
            self.lik.comps[modelP.comp0+i].model = modelP
            self.lik.comps[modelP.comp0+i].cweight = 0.0
            self.lik.comps[modelP.comp0+i].cmult = 0.0

    cdef void _decompModel(self, CxtLikModel *modelP) except *:
        # Normalize Pi as needed, and decompose Q.
        if CxLikQDecomp(self.lik.dim, modelP.rTri, modelP.piDiag, \
          modelP.piDiagNorm, &modelP.qNorm, modelP.qEigVecCube, \
          modelP.qEigVals):
            raise ValueError("Error decomposing Q")

        modelP.decomp = False

    cdef void _deallocModel(self, CxtLikModel *modelP, unsigned model):
        if modelP.rclass != NULL:
            free(modelP.rclass)
            modelP.rclass = NULL
        if modelP.rTri != NULL:
            free(modelP.rTri)
            modelP.rTri = NULL
        if modelP.piDiag != NULL:
            free(modelP.piDiag)
            modelP.piDiag = NULL
        if modelP.piDiagNorm != NULL:
            free(modelP.piDiagNorm)
            modelP.piDiagNorm = NULL
        if modelP.qEigVecCube != NULL:
            free(modelP.qEigVecCube)
            modelP.qEigVecCube = NULL
        if modelP.qEigVals != NULL:
            free(modelP.qEigVals)
            modelP.qEigVals = NULL

        if modelP.comp0 + modelP.clen < self.lik.compsLen:
            # Fill the hole in the comps vector.
            for model < i < self.lik.modelsLen:
                self.lik.models[i].comp0 -= modelP.clen
            memmove(&self.lik.comps[modelP.comp0], \
              &self.lik.comps[modelP.comp0+modelP.clen], \
              (self.lik.compsLen - (modelP.comp0+modelP.clen)) * \
              sizeof(CxtLikComp))

        self.lik.compsLen -= modelP.clen
        free(modelP)

    IF @enable_mpi@:
        cdef void configMpi(self, mpi.MPI_Comm mpiComm) except *:
            """
                Stripe data across the MPI nodes in mpiComm.
            """
            cdef int mpiSize, mpiRank
            cdef unsigned npad, nchars, mschars
            cdef double *stripeLnL, *siteLnL

            # Get communicator info.
            mpi.MPI_Comm_size(mpiComm, &mpiSize)
            mpi.MPI_Comm_rank(mpiComm, &mpiRank)
            # If necessary, modify padding to make the alignment width a
            # multiple of the MPI node count.  Note that multi-threaded
            # striping is disabled whenever MPI striping is enabled, because
            # there is no guarantee that each MPI node has the same number of
            # CPUs, thus making it impossible to base striping (and padding) on
            # the number of CPUs.
            if (self.alignment.nchars-self.alignment.npad) % mpiSize != 0:
                npad = mpiSize - \
                  ((self.alignment.nchars-self.alignment.npad) % mpiSize)
            else:
                npad = 0
            mschars = ((self.alignment.nchars-self.alignment.npad) + npad) / \
              mpiSize
            # Add padding to the alignment if necessary.
            if npad > self.alignment.npad:
                self.alignment.pad(self.char_.val2code(self.char_.any), \
                  npad - self.lik.npad)
                self.lik.charFreqs = self.alignment.freqs
            else:
                npad = self.lik.npad
            nchars = self.alignment.nchars
            # Resize siteLnL and stripeLnL before storing computed results.
            siteLnL = <double *>realloc(self.lik.siteLnL, nchars * \
              sizeof(double))
            if siteLnL == NULL:
                raise MemoryError("Error allocating siteLnL")
            self.lik.siteLnL = siteLnL
            stripeLnL = <double *>realloc(self.lik.stripeLnL, sizeof(double))
            if stripeLnL == NULL:
                raise MemoryError("Error allocating stripeLnL")
            self.lik.stripeLnL = stripeLnL
            # Store results.
            self.lik.mpiComm = mpiComm
            self.lik.mpiSize = mpiSize
            self.lik.mpiRank = mpiRank
            self.lik.nchars = nchars
            self.lik.npad = npad
            self.lik.mschars = mschars
            self.lik.cbase = mpiRank * self.lik.mschars
            self.lik.stripeWidth = mschars
            self.lik.nstripes = 1
            # Discard all CL data, since it is probably incorrectly sized.
            self.tree.clearAux()

        def enableMpi(self):
            """
                Stripe data across all MPI nodes.
            """
            self.configMpi(mpi.MPI_COMM_WORLD)

    cpdef Lik unpickle(self, str pickle):
        """
            Unpickle the partial Lik encoded by 'pickle', and fill in the
            missing details by copying them from 'self'.
        """
        cdef Lik ret

        ret = cPickle.loads(pickle)

        # Fill in missing details.
        ret._init2(self.alignment, self.char_)

        return ret

    cdef void _dup(self, Lik lik) except *:
        cdef unsigned i, m, ncat
        cdef CxtLikModel *toP, *frP

        assert lik.lik.modelsLen == 0

        for 0 <= i < self.lik.modelsLen:
            frP = self.lik.models[i]
            ncat = frP.clen
            if frP.invar:
                ncat -= 1
            m = lik.addModel(frP.weight, ncat, frP.catMedian, frP.invar)
            toP = lik.lik.models[m]
            assert toP.decomp
            assert toP.weight == frP.weight
            toP.rmult = frP.rmult
            memcpy(toP.rclass, frP.rclass, self.lik.rlen * sizeof(unsigned))
            memcpy(toP.rTri, frP.rTri, self.lik.rlen * sizeof(double))
            memcpy(toP.piDiag, frP.piDiag, self.lik.dim * sizeof(double))
            memcpy(toP.piDiagNorm, frP.piDiagNorm, self.lik.dim * \
              sizeof(double))

            toP.alpha = frP.alpha
            toP.catMedian = frP.catMedian

        for 0 <= i < self.lik.compsLen:
            lik.lik.comps[i].cweight = self.lik.comps[i].cweight
            lik.lik.comps[i].cmult = self.lik.comps[i].cmult

        IF @enable_mpi@:
            cdef double *siteLnL, *stripeLnL

            if self.lik.mpiComm != mpi.MPI_COMM_NULL:
                siteLnL = <double *>realloc(lik.lik.siteLnL, self.lik.nchars * \
                  sizeof(double))
                if siteLnL == NULL:
                    raise MemoryError("Error allocating siteLnL")
                lik.lik.siteLnL = siteLnL

                stripeLnL = <double *>realloc(lik.lik.stripeLnL, \
                  self.lik.mpiSize * sizeof(double))
                if stripeLnL == NULL:
                    raise MemoryError("Error allocating stripeLnL")
                lik.lik.stripeLnL = stripeLnL

                lik.lik.mpiComm = self.lik.mpiComm
                lik.lik.mpiSize = self.lik.mpiSize
                lik.lik.mpiRank = self.lik.mpiRank
                lik.lik.nchars = self.lik.nchars
                lik.lik.npad = self.lik.npad
                lik.lik.mschars = self.lik.mschars
                lik.lik.cbase = self.lik.cbase
                lik.lik.stripeWidth = self.lik.stripeWidth
                lik.lik.nstripes = self.lik.nstripes

    cpdef Lik dup(self):
        """
            Create a duplicate Lik object that can be manipulated independently
            of the original (with the exception that they refer to the same
            Alignment).  Copy all model parameters, and create a separate
            identical tree.  Do not copy cached conditional likelihood data;
            they can be recomputed if necessary.
        """
        cdef Lik ret
        cdef Tree tree

        ret = Lik()
        tree = self.tree.dup()
        ret._init0(tree)
        ret._init1(tree, self.lik.nchars, self.lik.dim, 0)
        ret._init2(self.alignment, self.char_)
        self._dup(ret)

        return ret

    cdef list _simulateRoot(self, list brks):
        cdef list seq, cdf
        cdef unsigned c, i, comp, brk
        cdef double pCum

        seq = []
        c = 0
        assert len(brks) == self.lik.compsLen
        for 0 <= comp < self.lik.compsLen:
            brk = brks[comp]
            # Create a CDF for state frequencies.
            cdf = []
            pCum = 0.0
            for 0 <= i < self.lik.dim-1:
                pCum += self.lik.comps[comp].model.piDiagNorm[i]
                if pCum > 1.0:
                    pCum = 1.0
                cdf.append(pCum)
            cdf.append(1.0)
            # Generate random characters, proportional to the state frequencies
            # for this model component.
            while c < brk:
                u = random.random()
                for 0 <= i < self.lik.dim:
                    if cdf[i] >= u:
                        seq.append(i)
                        break
                assert i < self.lik.dim
                c += 1
                assert len(seq) == c

        return seq

    cdef list _simulateChild(self, list brks, list parSeq, double brlen):
        cdef list seq, cdfs, cdf
        cdef double *P
        cdef unsigned c, i, j, comp, brk
        cdef double pCum

        P = <double *>malloc(self.lik.dim * self.lik.dim * sizeof(double))
        if P == NULL:
            raise MemoryError("Error allocating P")
        try:
            seq = []
            c = 0
            assert len(brks) == self.lik.compsLen
            for 0 <= comp < self.lik.compsLen:
                brk = brks[comp]
                # Create a list of CDFs for state transition probabilities,
                # one for each starting state.
                cdfs = []
                CxLikPt(self.lik.dim, P, \
                  self.lik.comps[comp].model.qEigVecCube,
                  self.lik.comps[comp].model.qEigVals, brlen)
                for 0 <= i < self.lik.dim:
                    cdf = []
                    pCum = 0.0
                    for 0 <= j < self.lik.dim-1:
                        pCum += P[i*self.lik.dim + j]
                        if pCum > 1.0:
                            pCum = 1.0
                        cdf.append(pCum)
                    cdf.append(1.0)
                    cdfs.append(cdf)
                # Generate random characters, based on starting state and
                # transition probabilities.
                while c < brk:
                    u = random.random()
                    cdf = cdfs[parSeq[c]]
                    for 0 <= i < self.lik.dim:
                        if cdf[i] >= u:
                            seq.append(i)
                            break
                    assert i < self.lik.dim
                    c += 1
                    assert len(seq) == c
        finally:
            free(P)

        return seq

    cdef void _simulateRecurse(self, str i2c, list brks, list parSeq, \
      Ring ring) except *:
        cdef list seq
        cdef Taxon taxon
        cdef Ring r

        seq = self._simulateChild(brks, parSeq, ring.edge.length)
        taxon = ring.node.getTaxon()
        if taxon is not None:
            self.alignment.setSeq(self.alignment.taxaMap.indGet(taxon), 0, \
              "".join([i2c[seq[i]] for i in xrange(len(seq))]))
        for r in ring.siblings():
            self._simulateRecurse(i2c, brks, seq, r.other)

    cdef void _simulate(self) except *:
        cdef list brks, seq
        cdef double wSum
        cdef unsigned i
        cdef Node base
        cdef Taxon taxon
        cdef Ring ring, r
        cdef str i2c

        # Create a lookup table from state index to character.  For DNA this is
        # "ACGT".
        i2c = "".join(self.char_.pcodes())

        # Compute the breaks between characters simulated under each model
        # component.  For example, given a 1000-character alignment and the
        # following component weights, breaks will look like:
        #
        #   weightScaled:  0.25, 0.10, 0.00, 0.00, 0.60, 0.05
        #   brks        : [   0,  250,  250,  250,  950, 1000]
        self.prep()
        brks = []
        wSum = 0.0
        for 0 <= i < self.lik.compsLen-1:
            wSum += self.lik.comps[i].weightScaled
            if wSum > 1.0:
                wSum = 1.0
            brks.append(<int>ceil(wSum * \
              <double>(self.alignment.nchars - self.alignment.npad)))
        brks.append(self.alignment.nchars - self.alignment.npad)

        seq = self._simulateRoot(brks)
        base = self.tree.getBase()
        taxon = base.getTaxon()
        if taxon is not None:
            self.alignment.setSeq(self.alignment.taxaMap.indGet(taxon), 0, \
              "".join([i2c[seq[i]] for i in xrange(len(seq))]))
        ring = base.ring
        if ring is not None:
            for r in ring:
                self._simulateRecurse(i2c, brks, seq, r.other)

        # Clean up.
        self.alignment.deCol()

    cpdef Lik simulate(self, unsigned nchars=0):
        """
            Create a derivative Lik instance that is a duplicate (see the dup()
            method), except that its alignment is replaced with one that is
            filled with simulated sequences, based on the current model
            parameters (including tree topology and branch lengths).  If nchars
            is unspecified, the simulated alignment has the same number of
            characters as the original alignment.
        """
        cdef Lik ret
        cdef Tree Tree
        cdef unsigned npad
        cdef Alignment alignment

        ret = Lik()
        tree = self.tree.dup()
        if nchars == 0:
            nchars = self.alignment.nchars - self.alignment.npad
        alignment = Alignment(None, None, self.alignment.taxaMap, nchars, \
          self.alignment.charType, True, True)
        # Pad the alignment if its width isn't a multiple of the stripe
        # width.
        npad = self._computeNpad(nchars, self._computeStripeWidth(nchars))
        if npad != 0:
            alignment.pad(self.char_.val2code(self.char_.any), npad)
            nchars += npad
        ret._init0(tree)
        ret._init1(tree, nchars, self.lik.dim, 0)
        ret._init2(alignment, self.char_)
        self._dup(ret)

        ret._simulate()
        return ret

    cpdef Lik clone(self):
        """
            Return a Lik instance that has exactly the same model parameters,
            and shares the same underlying tree.  Cached likelihood data are
            not copied, since the intent of clone() is to provide a separate
            Lik in which to modify model parameters such that the cached data
            would be invalidated anyway.

            Internally, each Lik has at most one mated Lik instance that is
            used for clone().  This means that in the following code, lik0 and
            lik2 are the same instance:

              lik0 = Lik(...)
              lik1 = lik0.clone()
              lik2 = lik1.clone()
        """
        cdef Lik ret
        cdef int i
        cdef CxtLikModel *modelP
        cdef bint resize

        if self.mate is not None:
            ret = self.mate

            ret.lik.invalidate = True
            # There's no need to discard internal-node cLMat's unless compsLen
            # differs between mates.
            resize = (ret.lik.compsLen != self.lik.compsLen)

            # Clear ret's models/comps vectors.  Iterate downward to avoid
            # gratuitous memory moves.
            for ret.lik.modelsLen > i >= 0:
                ret.delModel(i)
            assert ret.lik.modelsLen == 0
            assert ret.lik.compsLen == 0
        else:
            assert self.lik.polarity == 0
            resize = False
            ret = Lik()
            ret._init1(self.tree, self.lik.nchars, self.lik.dim, 1)
            ret._init2(self.alignment, self.char_)
            ret.mate = self
            self.mate = ret

        self._dup(ret)
        ret.lik.resize = resize
        ret.lik.reweight = True

        return ret

    cpdef double getWNorm(self) except -1.0:
        """
            Get the weighted normalization factor that is used to normalize all
            model Q matrices, so that the mixture's mean relative mutation
            rate is 1.0.
        """
        return self.lik.wNorm

    cpdef unsigned nmodels(self):
        """
            Return the current number of models being mixed.  Models can be
            referenced by their offset within the mixture model vector.
        """
        return self.lik.modelsLen

    cpdef unsigned addModel(self, double weight, unsigned ncat=1, \
      bint catMedian=False, bint invar=False) except *:
        """
            Append a model to the mixture, and return its index in the model
            vector.
        """
        cdef unsigned ret, i
        cdef CxtLikModel **models

        assert ncat > 0

        if self.lik.modelsLen == self.lik.modelsMax:
            models = <CxtLikModel **>realloc(self.lik.models, \
              (self.lik.modelsMax + 1) * sizeof(CxtLikModel *))
            if models == NULL:
                raise MemoryError("Error reallocating models")
            self.lik.models = models
            self.lik.modelsMax += 1
        self.lik.models[self.lik.modelsLen] = self._allocModel(ncat, invar)
        self._initModel(self.lik.models[self.lik.modelsLen], weight, \
          catMedian, invar)
        ret = self.lik.modelsLen
        self.lik.modelsLen += 1

        self.lik.resize = True

        return ret

    cpdef delModel(self, unsigned model):
        """
            Remove the model from the mixture vector, and renumber remaining
            models so that model numbers remain contiguous.
        """
        cdef CxtLikModel *modelP

        assert model < self.lik.modelsLen
        modelP = self.lik.models[model]

        if self.lik.modelsLen == 0:
            raise ValueError("No models remain")

        self._deallocModel(modelP, model)

        self.lik.resize = True

        if model != self.lik.modelsLen - 1:
            # Fill the hole in the models vector.
            memmove(&self.lik.models[model], &self.lik.models[model+1], \
              (self.lik.modelsLen - (model+1)) * sizeof(CxtLikModel *))
        self.lik.modelsLen -= 1

    cpdef double getWeight(self, unsigned model) except -1.0:
        """
            Get the relative weight for the specified model in the mixture
            vector.
        """
        cdef CxtLikModel *modelP

        assert model < self.lik.modelsLen
        modelP = self.lik.models[model]

        return modelP.weight

    cpdef setWeight(self, unsigned model, double weight):
        """
            Set the relative weight for the specified model in the mixture
            vector.
        """
        cdef CxtLikModel *modelP

        assert model < self.lik.modelsLen
        modelP = self.lik.models[model]
        assert weight >= 0.0

        modelP.weight = weight
        self.lik.reweight = True

    cpdef double getRmult(self, unsigned model) except -1.0:
        """
            Get the rate multiplier for the specified model.
        """
        cdef CxtLikModel *modelP

        assert model < self.lik.modelsLen
        modelP = self.lik.models[model]

        return modelP.rmult

    cpdef setRmult(self, unsigned model, double rmult):
        """
            Set the rate multiplier for the specified model.
        """
        cdef CxtLikModel *modelP

        assert model < self.lik.modelsLen
        modelP = self.lik.models[model]
        assert rmult >= 0.0

        modelP.rmult = rmult
        self.lik.reweight = True

    cpdef list getRclass(self, unsigned model):
        """
            Return a list of integers that encodes relative mutation rate class
            partitioning for the specified model.  For DNA, where the R matrix
            is
               _            _
              |  -  a  b  c  |
              |              |
              |  a  -  d  e  |
              |              | ,
              |  b  d  -  f  |
              |              |
              |_ c  e  f  - _|

            the return list is [a,b,c,d,e,f].  GTR would be [0,1,2,3,4,5].
        """
        cdef CxtLikModel *modelP
        cdef unsigned i

        assert model < self.lik.modelsLen
        modelP = self.lik.models[model]

        return [modelP.rclass[i] for i in xrange(self.lik.rlen)]

    cpdef setRclass(self, unsigned model, list rclass, list rates=None):
        """
            Set the relative mutation rate class partitioning for the specified
            model.  For DNA, where the R matrix is
               _            _
              |  -  a  b  c  |
              |              |
              |  a  -  d  e  |
              |              | ,
              |  b  d  -  f  |
              |              |
              |_ c  e  f  - _|

            rclass=[0,1,2,3,4,5] would set the model to GTR.  Note that rate
            class integer lists must be in a canonized form, such that the first
            element of the list is always 0, and subsequent elements are at most
            one greater than all preceeding elements.

            The relative mutation rates can be specified via the 'rates'
            parameter, which if omitted causes all rates to be set to 1.

            Named model rclass equivalents:

              base freqs   |
              ==    !=     |
              .............|
              model        | rclass
              -------------+--------------
              JC     F81   | [0,0,0,0,0,0]
              K80    HKY   | [0,1,0,0,1,0]
              TrNef  TrN   | [0,1,0,0,2,0]
              K81    K81uf | [0,1,2,2,1,0]
              TIMef  TIM   | [0,1,2,2,3,0]
                     TVM   | [0,1,2,3,1,4]
              SYM    GTR   | [0,1,2,3,4,5]
        """
        cdef CxtLikModel *modelP
        cdef unsigned rMax, r, i

        assert model < self.lik.modelsLen
        modelP = self.lik.models[model]

        rMax = 0
        for 0 <= i < self.lik.rlen:
            r = rclass[i]
            if r == rMax + 1:
                rMax += 1
            elif r > rMax + 1:
                raise ValueError("Invalid rclass specification")

        for 0 <= i < self.lik.rlen:
            modelP.rclass[i] = rclass[i]

        if rates is None:
            for 0 <= i < self.lik.rlen:
                modelP.rTri[i] = 1.0
        else:
            if len(rates) != rMax + 1:
                raise ValueError("Incorrect rates list length")
            for 0 <= i < self.lik.rlen:
                if type(rates[modelP.rclass[i]]) is not float:
                    raise ValueError("Incorrect type for rate %d (%r)" % \
                      (rates[modelP.rclass[i]], i))
            for 0 <= i < self.lik.rlen:
                modelP.rTri[i] = rates[modelP.rclass[i]]

        modelP.decomp = True

    cpdef unsigned getNrates(self, unsigned model) except 0:
        """
            Get the number of relative mutation rate classes for the specified
            model.
        """
        cdef unsigned nrates, i
        cdef CxtLikModel *modelP

        assert model < self.lik.modelsLen
        modelP = self.lik.models[model]

        nrates = 1
        for 0 <= i < self.lik.rlen:
            rInd = modelP.rclass[i]
            if rInd + 1 > nrates:
                nrates = rInd + 1

        return nrates

    cpdef double getRate(self, unsigned model, unsigned i) except -1.0:
        """
            Get the value of the i'th relative mutation rate class for the
            specified model.
        """
        cdef CxtLikModel *modelP
        cdef unsigned j

        assert model < self.lik.modelsLen
        modelP = self.lik.models[model]
        assert i < self.lik.rlen

        for 0 <= j < self.lik.rlen:
            if modelP.rclass[j] == i:
                return modelP.rTri[j]

        raise ValueError("Invalid rate class index")

    cpdef setRate(self, unsigned model, unsigned i, double rate):
        """
            Set the value of the i'th relative mutation rate class for the
            specified model.
        """
        cdef CxtLikModel *modelP
        cdef unsigned j
        cdef bint iValid, decomp

        assert model < self.lik.modelsLen
        modelP = self.lik.models[model]
        assert i < self.lik.rlen
        assert rate >= 0.0

        iValid = False
        decomp = False
        for 0 <= j < self.lik.rlen:
            if modelP.rclass[j] == i:
                if modelP.rTri[j] != rate:
                    decomp = True
                modelP.rTri[j] = rate
                iValid = True
        if not iValid:
            raise ValueError("Invalid rate class index")
        if decomp:
            modelP.decomp = True

    cpdef double getFreq(self, unsigned model, unsigned i) except -1.0:
        """
            Get the frequency parameter for state i, for the specified model.
        """
        cdef CxtLikModel *modelP

        assert model < self.lik.modelsLen
        modelP = self.lik.models[model]
        assert i < self.lik.dim

        return modelP.piDiag[i]

    cpdef setFreq(self, unsigned model, unsigned i, double freq):
        """
            Set the frequency parameter for state i, for the specified model.
        """
        cdef CxtLikModel *modelP

        assert model < self.lik.modelsLen
        modelP = self.lik.models[model]
        assert i < self.lik.dim
        assert freq > 0.0

        if modelP.piDiag[i] != freq:
            modelP.piDiag[i] = freq
            modelP.decomp = True

    cpdef double getAlpha(self, unsigned model):
        """
            Get the alpha shape parameter controlling Gamma-distributed
            mutation rate variation, for the specified model.
        """
        cdef CxtLikModel *modelP

        assert model < self.lik.modelsLen
        modelP = self.lik.models[model]

        return modelP.alpha

    cpdef setAlpha(self, unsigned model, double alpha):
        """
            Set the alpha shape parameter controlling Gamma-distributed
            mutation rate variation, for the specified model.
        """
        cdef CxtLikModel *modelP
        cdef double lnGammaA, lnGammaA1, pt, sum, cweight
        cdef unsigned ncat, i

        assert model < self.lik.modelsLen
        modelP = self.lik.models[model]
        ncat = modelP.clen
        if modelP.invar:
            ncat -= 1
        assert ncat > 1

        if modelP.alpha != alpha:
            modelP.alpha = alpha
            if alpha != INFINITY:
                # Discretize Gamma-distributed rates.
                lnGammaA = lnGamma(alpha)
                if not modelP.catMedian: # Category means.
                    # Use properly scaled ptChi2() results to compute boundaries
                    # between equal-size rate categories.  Then use the
                    # boundaries to compute category means.
                    lnGammaA1 = lnGamma(alpha+1.0)
                    for 0 <= i < ncat - 1:
                        pt = ptChi2((<double>i+1.0) / <double>ncat, alpha*2.0, \
                          lnGammaA)
                        if pt == -1.0:
                            raise OverflowError(\
                              "Error discretizing gamma (ncat=%d, alpha=%e)" % \
                              (ncat, alpha))
                        self.lik.comps[modelP.comp0+i].cmult = \
                          gammaI(pt/2.0, alpha+1.0, lnGammaA1)
                    self.lik.comps[modelP.comp0+ncat-1].cmult = 1.0

                    # Convert to relative rates and rescale to a mean rate of 1.
                    for ncat - 1 >= i > 0:
                        self.lik.comps[modelP.comp0+i].cmult -= \
                          self.lik.comps[modelP.comp0+i-1].cmult
                        self.lik.comps[modelP.comp0+i].cmult *= <double>ncat
                    self.lik.comps[modelP.comp0].cmult *= <double>ncat
                else: # Category medians.
                    # Use properly scaled ptChi2() results to compute point
                    # values for the medians of ncat equal-sized partitions,
                    sum = 0.0
                    for 0 <= i < ncat:
                        pt = ptChi2((<double>(i*2)+1.0) / \
                          <double>(ncat*2), alpha*2.0, lnGammaA)
                        if pt == -1.0:
                            raise OverflowError(\
                              "Error discretizing gamma (ncat=%d, alpha=%e)" % \
                              (ncat, alpha))
                        pt /= alpha * 2.0
                        sum += pt
                        self.lik.comps[modelP.comp0+i].cmult = pt

                    # Rescale to a mean rate of 1.
                    for 0 <= i < ncat:
                        self.lik.comps[modelP.comp0+i].cmult *= \
                          <double>ncat / sum

                # Evenly weight the +G components.
                cweight = 0.0
                for 0 <= i < ncat:
                    cweight += self.lik.comps[modelP.comp0+i].cweight
                cweight /= <double>ncat
                for 0 <= i < ncat:
                    self.lik.comps[modelP.comp0+i].cweight = cweight
            else:
                # No Gamma-distributed rates.
                # Use only the first +G component.
                cweight = 0.0
                for 0 <= i < ncat:
                    cweight += self.lik.comps[modelP.comp0+i].cweight
                self.lik.comps[modelP.comp0].cweight = cweight
                self.lik.comps[modelP.comp0].cmult = 1.0
                for 1 <= i < ncat:
                    self.lik.comps[modelP.comp0+i].cweight = 0.0
                    self.lik.comps[modelP.comp0+i].cmult = 0.0

            self.lik.invalidate = True

    cpdef unsigned getNcat(self, unsigned model) except *:
        """
            Get the number of discrete rates for Gamma-distributed mutation
            rate variation, for the specified model.
        """
        cdef unsigned ncat
        cdef CxtLikModel *modelP

        assert model < self.lik.modelsLen
        modelP = self.lik.models[model]

        ncat = modelP.clen
        if modelP.invar:
            ncat -= 1
        return ncat

    cpdef bint getCatMedian(self, unsigned model) except *:
        """
            Get the Gamma-distributed rate discretization algorithm for the
            specified model.  If False, category means are used; if True,
            category medians are used.
        """
        cdef CxtLikModel *modelP

        assert model < self.lik.modelsLen
        modelP = self.lik.models[model]

        return modelP.catMedian

    cpdef bint getInvar(self, unsigned model) except *:
        """
            Get whether the specified model supports a proportion of invariable
            sites (i.e. whether it is a +I model).
        """
        cdef CxtLikModel *modelP

        assert model < self.lik.modelsLen
        modelP = self.lik.models[model]

        return modelP.invar

    cpdef double getWVar(self, unsigned model) except *:
        """
            Get the relative weight for variable sites.
        """
        cdef double wVar
        cdef CxtLikModel *modelP
        cdef unsigned ncat, i

        assert model < self.lik.modelsLen
        modelP = self.lik.models[model]

        if modelP.alpha == INFINITY:
            wVar = self.lik.comps[modelP.comp0].cweight
        else:
            ncat = modelP.clen
            if modelP.invar:
                ncat -= 1

            wVar = 0.0
            for 0 <= i < ncat:
                wVar += self.lik.comps[modelP.comp0+i].cweight

        return wVar

    cpdef setWVar(self, unsigned model, double wVar):
        """
            Set the relative weight for variable sites.
        """
        cdef CxtLikModel *modelP
        cdef unsigned ncat, i

        assert model < self.lik.modelsLen
        modelP = self.lik.models[model]
        assert wVar >= 0.0

        if modelP.alpha == INFINITY:
            self.lik.comps[modelP.comp0].cweight = wVar
        else:
            ncat = modelP.clen
            if modelP.invar:
                ncat -= 1

            for 0 <= i < ncat:
                self.lik.comps[modelP.comp0+i].cweight = wVar / <double>ncat
        self.lik.reweight = True

    cpdef double getWInvar(self, unsigned model) except *:
        """
            Get the relative weight for invariable sites.
        """
        cdef CxtLikModel *modelP

        assert model < self.lik.modelsLen
        modelP = self.lik.models[model]
        if not modelP.invar:
            return 0.0
        else:
            return self.lik.comps[modelP.comp0+modelP.clen-1].cweight

    cpdef setWInvar(self, unsigned model, double wInvar):
        """
            Set the relative weight for invariable sites.
        """
        cdef CxtLikModel *modelP

        assert model < self.lik.modelsLen
        modelP = self.lik.models[model]
        if not modelP.invar and wInvar != 0.0:
            raise ValueError("Model does not support invariable sites.")
        assert wInvar >= 0.0

        self.lik.comps[modelP.comp0+modelP.clen-1].cweight = wInvar
        self.lik.reweight = True

    cdef void _planAppend(self, CxeLikStep variant, unsigned ntrail, \
      CL parentCL, CL childCL, double edgeLen) except *:
        cdef CxtLikStep *step

        assert parentCL is not None
        assert childCL is not None
        assert ntrail != 0 or variant in (CxeLikStepMergeL, CxeLikStepMergeI)

        step = &self.lik.steps[self.lik.stepsLen]
        self.lik.stepsLen += 1

        step.variant = variant
        step.ntrail = ntrail
        step.parentCL = &parentCL.cLs[self.lik.polarity]
        assert step.parentCL != NULL
        assert step.parentCL.cLMat != NULL
        assert step.parentCL.lnScale != NULL
        if childCL.cLs[1].cLMat == NULL:
            # Be careful with leaf nodes to always use the first (and only)
            # cLs element.
            assert childCL.cLs[0].cLMat != NULL
            step.childCL = &childCL.cLs[0]
        else:
            assert childCL.cLs[self.lik.polarity].cLMat != NULL
            step.childCL = &childCL.cLs[self.lik.polarity]
        assert step.childCL != NULL
        assert step.childCL.cLMat != NULL
        assert step.childCL.lnScale != NULL
        if edgeLen < 0.0:
            raise ValueError("Negative branch length")
        step.edgeLen = edgeLen

    cdef void _planRecurse(self, Ring ring, CL parent, unsigned nSibs,
      double edgeLen) except *:
        cdef CL cL, flushCL
        cdef Taxon taxon
        cdef unsigned degree, i, j
        cdef char *chars
        cdef double *cLMat
        cdef int ind, val
        cdef Ring r
        cdef CxtLikCL *cLC, *pCLC

        assert edgeLen == ring.edge.length or edgeLen == 0.0

        cL = <CL>ring.aux
        degree = ring.node.getDegree()

        # Initialize/expand cL as necessary.  Leaf nodes are initialized here
        # rather than in the Lik constructor so that sequential addition
        # methods will work correctly.
        if degree == 1:
            if cL is None:
                taxon = ring.node.taxon
                if taxon is None:
                    raise ValueError("Leaf node missing taxon")
                # Leaf nodes only need cLs[0], since character data can be
                # shared by all model components.
                cL = CL()
                cL.prepare(0, self.lik.mschars, self.lik.dim, 1)
                # Set lnScale entries to 0.0 for leaves.  This is the only
                # place that explicit initialization of lnScale is necessary.
                memset(cL.cLs[0].lnScale, 0, self.lik.mschars * sizeof(double))
                ring.aux = cL

                ind = self.alignment.taxaMap.indGet(taxon)
                if ind == -1:
                    raise ValueError( \
                      "Taxon %r missing from alignment's taxa map" % \
                      taxon.label)
                chars = self.alignment.getRow(ind)
                cLMat = cL.cLs[0].cLMat
                for 0 <= i < self.lik.mschars:
                    val = self.char_.code2val(chr(chars[self.lik.cbase + i]))
                    if val == 0:
                        val = self.char_.any
                    for 0 <= j < self.lik.dim:
                        if val & (1 << j):
                            cLMat[i*self.lik.dim + j] = 1.0
                        else:
                            cLMat[i*self.lik.dim + j] = 0.0
        else:
            if cL is None:
                cL = CL()
                cL.prepare(self.lik.polarity, self.lik.mschars, self.lik.dim, \
                  self.lik.compsLen)
                ring.aux = cL
            else:
                if self.lik.resize:
                    cL.resize(self.lik.polarity, self.lik.mschars, \
                      self.lik.dim, self.lik.compsLen)
                else:
                    cL.prepare(self.lik.polarity, self.lik.mschars, \
                      self.lik.dim, self.lik.compsLen)

        # Recurse.
        for r in ring.siblings():
            if self.lik.invalidate:
                # Try to clean up invalid CL caches.
                flushCL = <CL>r.aux
                if flushCL is not None:
                    flushCL.flush(self.lik.polarity)

            self._planRecurse(r.other, cL, degree, r.edge.length)

        # Check whether the current tree topology is compatible with the
        # parent's cache.  If not, invalidate the parent's cache.
        cLC = &cL.cLs[self.lik.polarity]
        pCLC = &parent.cLs[self.lik.polarity]
        if self.lik.invalidate or cLC.parent is not pCLC or cLC.nSibs != nSibs \
          or cLC.edgeLen != edgeLen:
            pCLC.valid = False
            cLC.parent = &parent.cLs[self.lik.polarity]
            cLC.nSibs = nSibs
            cLC.edgeLen = edgeLen

        # Check cache validity.
        if degree != 1:
            if self.lik.invalidate or (not cLC.valid):
                r = ring.next
                for 1 <= j < degree:
                    assert r != ring
                    if r.other.node.getDegree() > 1:
                        if j == 1:
                            variant = CxeLikStepComputeI
                        else:
                            variant = CxeLikStepMergeI
                    else:
                        if j == 1:
                            variant = CxeLikStepComputeL
                        else:
                            variant = CxeLikStepMergeL
                    self._planAppend(variant, degree-j-1, cL, <CL>r.other.aux, \
                      r.edge.length)
                    r = r.next
                cLC.valid = True
                # Propagate invalidation to the parent.
                pCLC.valid = False

    cdef void _plan(self, Node root) except *:
        cdef Ring ring, r
        cdef unsigned degree, ntrail
        cdef CL flushCL
        cdef CxeLikStep variant

        if self.lik.resize:
            self.rootCL.resize(self.lik.polarity, self.lik.mschars, \
              self.lik.dim, self.lik.compsLen)
        else:
            self.rootCL.prepare(self.lik.polarity, self.lik.mschars, \
              self.lik.dim, self.lik.compsLen)

        if root is None:
            root = self.tree.base
        degree = root.getDegree()
        ring = root.ring
        if degree == 1:
            # The root is a leaf node, which requires special handling.  Treat
            # the root as if it were a child separated from the root by a
            # 0-length branch.
            self._planRecurse(ring.other, self.rootCL, 2, ring.edge.length)
            self._planRecurse(ring, self.rootCL, 2, 0.0)

            if self.lik.invalidate or \
              (not self.rootCL.cLs[self.lik.polarity].valid):
                if ring.other.node.getDegree() > 1:
                    variant = CxeLikStepComputeI
                else:
                    variant = CxeLikStepComputeL
                self._planAppend(variant, 1, self.rootCL, <CL>ring.other.aux, \
                  ring.edge.length)
                self._planAppend(CxeLikStepMergeL, 0, self.rootCL, \
                  <CL>ring.aux, 0.0)
                self.rootCL.cLs[self.lik.polarity].valid = True
        else:
            for r in ring:
                if self.lik.invalidate:
                    # Try to clean up invalid CL caches.
                    flushCL = <CL>r.aux
                    if flushCL is not None:
                        flushCL.flush(self.lik.polarity)

                self._planRecurse(r.other, self.rootCL, degree, r.edge.length)

            if self.lik.invalidate or \
              (not self.rootCL.cLs[self.lik.polarity].valid):
                if ring.other.node.getDegree() > 1:
                    variant = CxeLikStepComputeI
                else:
                    variant = CxeLikStepComputeL
                ntrail = degree - 1
                self._planAppend(variant, ntrail, self.rootCL, \
                  <CL>ring.other.aux, ring.edge.length)
                for r in ring.siblings():
                    if r.other.node.getDegree() > 1:
                        variant = CxeLikStepMergeI
                    else:
                        variant = CxeLikStepMergeL
                    ntrail -= 1
                    self._planAppend(variant, ntrail, self.rootCL, \
                      <CL>r.other.aux, r.edge.length)
                self.rootCL.cLs[self.lik.polarity].valid = True

        # Now that execution planning is complete, clear flags that have been
        # acted on.
        self.lik.resize = False
        self.lik.invalidate = False

    cpdef prep(self):
        """
            Recompute all internal model parameters necessary for lnL
            computation.
        """
        cdef double wSum, wSumC, qMean, wNorm
        cdef unsigned i, j
        cdef CxtLikModel *modelP
        cdef CxtLikComp *comp

        if self.lik.modelsLen == 0:
            raise ValueError("Empty model mixture")

        if self.lik.resize:
            self.lik.reweight = True

        if self.lik.reweight:
            # Scale weights.
            wSum = 0.0
            for 0 <= i < self.lik.modelsLen:
                modelP = self.lik.models[i]
                wSum += modelP.weight
            if wSum == 0.0:
                raise ValueError("At least one model must have non-zero weight")

            for 0 <= i < self.lik.modelsLen:
                modelP = self.lik.models[i]
                wSumC = 0.0
                for 0 <= j < modelP.clen:
                    comp = &self.lik.comps[modelP.comp0+j]
                    wSumC += comp.cweight
                assert wSumC > 0.0
                for 0 <= j < modelP.clen:
                    comp = &self.lik.comps[modelP.comp0+j]
                    comp.weightScaled = (modelP.weight / wSum) * \
                      (comp.cweight / wSumC)
            self.lik.invalidate = True
            self.lik.reweight = False

        # Make sure that the mixture models are up to date.
        for 0 <= i < self.lik.modelsLen:
            modelP = self.lik.models[i]
            # Update model parameters, if necessary.
            if modelP.decomp:
                self._decompModel(modelP)

        # Recompute wNorm.
        qMean = 0.0
        for 0 <= i < self.lik.compsLen:
            comp = &self.lik.comps[i]
            if comp.weightScaled != 0.0:
                # Take special care to skip over +I components, since we don't
                # actually create Q matrices for them.  Were the Q matrices to
                # exist, qNorm would be INFINITY for them, and they would
                # contribute 0 to qMean.
                if comp.cmult != 0.0:
                    modelP = comp.model
                    qMean += comp.weightScaled * (modelP.rmult / modelP.qNorm)
        wNorm = 1.0 / qMean
        if wNorm != self.lik.wNorm:
            # Invalidate the cache, since wNorm impacts all CL's.
            self.lik.invalidate = True
            self.lik.wNorm = wNorm

    cdef void _prep(self, Node root) except *:
        cdef unsigned stepsMax
        cdef CxtLikStep *steps

        self.prep()

        # Expand steps, if necessary.
        stepsMax = (2 * self.alignment.ntaxa) - 2
        if self.lik.stepsMax < stepsMax:
            steps = <CxtLikStep *>realloc(self.lik.steps, stepsMax * \
              sizeof(CxtLikStep))
            if self.lik.steps == NULL:
                raise MemoryError("Error reallocating steps")
            self.lik.steps = steps
            self.lik.stepsMax = stepsMax
        self.lik.stepsLen = 0

        # Generate the execution plan via post-order tree traversal.
        self._plan(root)

    cpdef double lnL(self, Node root=None) except 1.0:
        """
            Compute the log-likelihood.  Use the tree's base node as the root
            for computation, unless a root is specified.
        """
        cdef double ret
        cdef unsigned i

        # Prepare data structures and compute the execution plan.
        self._prep(root)

        # Execute the plan.
        CxLikExecute(self.lik)

        # Sum stripe log-likelihoods.
        ret = 0.0
        IF @enable_mpi@:
            if self.lik.mpiComm != mpi.MPI_COMM_NULL:
                mpi.MPI_Allreduce(&self.lik.stripeLnL[0], &ret, 1, \
                  mpi.MPI_DOUBLE, mpi.MPI_SUM, self.lik.mpiComm)
            else:
                for 0 <= i < self.lik.nstripes:
                    ret += self.lik.stripeLnL[i]
        ELSE:
            for 0 <= i < self.lik.nstripes:
                ret += self.lik.stripeLnL[i]

        IF LikDebug:
            # Validate with a fresh Lik, in order to detect cache-related
            # flaws.
            cdef Lik lik = self.dup()
            lik._prep(root)
            CxLikExecute(lik.lik)
            cdef double lnL2 = 0.0
            IF @enable_mpi@:
                if lik.lik.mpiComm != mpi.MPI_COMM_NULL:
                    mpi.MPI_Allreduce(&lik.lik.stripeLnL[0], &lnL2, 1, \
                      mpi.MPI_DOUBLE, mpi.MPI_SUM, lik.lik.mpiComm)
                else:
                    for 0 <= i < lik.lik.nstripes:
                        lnL2 += lik.lik.stripeLnL[i]
            ELSE:
                for 0 <= i < lik.lik.nstripes:
                    lnL2 += lik.lik.stripeLnL[i]
            if not (0.99 < lnL2/ret and lnL2/ret < 1.01) and \
              not (isinf(lnL2) == -1 and isinf(ret) == -1):
                raise AssertionError( \
                  "Incorrect lnL: %f (should be approximately %f)" % \
                  (ret, lnL2))

        return ret

    cpdef list siteLnLs(self, Node root=None):
        """
            Compute the site log-likelihoods.  Use the tree's base node as the
            root for computation, unless a root is specified.
        """
        cdef list ret
        cdef double *lnLs

        # Prepare data structures and compute the execution plan.
        self._prep(root)

        # Execute the plan.
        CxLikExecute(self.lik)

        IF @enable_mpi@:
            if self.lik.mpiComm != mpi.MPI_COMM_NULL:
                mpi.MPI_Allgather(mpi.MPI_IN_PLACE, 0, mpi.MPI_DATATYPE_NULL, \
                  self.lik.siteLnL, self.lik.mschars, mpi.MPI_DOUBLE, \
                  self.lik.mpiComm)

        # Copy lnLs C array values into a Python list.
        ret = [self.lik.siteLnL[i] \
          for i in xrange(self.lik.nchars-self.lik.npad)]

        return ret

    cpdef flush(self):
        """
            Flush all cached conditional likelihood data.
        """
        cdef list edges
        cdef unsigned i
        cdef Edge edge
        cdef Ring ring
        cdef CL cL

        edges = self.tree.getEdges()
        for 0 <= i < len(edges):
            edge = <Edge>edges[i]
            ring = edge.ring
            if ring.node.getDegree() > 1:
                cL = <CL>ring.aux
                if cL is not None:
                    cL.flush(self.lik.polarity)
            ring = ring.other
            if ring.node.getDegree() > 1:
                cL = <CL>ring.aux
                if cL is not None:
                    cL.flush(self.lik.polarity)
