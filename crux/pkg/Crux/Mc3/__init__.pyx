"""
    Reversible jump Metropolis-coupled Markov chain Monte Carlo.

    Mc3 uses reversible jump Metropolis-coupled Markov chain Monte Carlo to
    create a representative sample (assuming convergence) of the stationary
    posterior distribution, based on the GTR+Gamma family of sequence evolution
    models.

    Polytomous tree topologies are sampled using the methods described by Lewis
    et al. (2005), and restrictions of the fully parameterized GTR model are
    sampled using the methods described by Huelsenbeck et al.  (2004).

    Extending tree bisection and reconnection (eTBR) is used to change topology
    and to modify branch lengths (Lakner et al. 2008).  This implementation is
    generalized to support polytomies, as well as to apply to all edges
    (including leaf edges).  The latter generalization allows topology
    tranformations even on trees that have a single internal edge.

    Models with/without Gamma-distributed relative mutation rate variation (+G)
    are sampled using novel methods.

    Mixtures (weighted Q matrices) of varying degree are sampled using novel
    methods.

    If MPI (http://www.mpi-forum.org/) support is enabled, the methods
    described by Altekar et al. (2004) are used to run chains in parallel.
    Ideally, the total number of chains (number of independent runs times
    number of Metropolis-coupled chains per run) should be an even multiple of
    the number of MPI nodes.

    Convergence (based on log-likelihoods) is monitored using the
    interval-based coverage ratio diagnostic described at the bottom of page
    441 in Brooks and Gelman (1998).  The diagnostic takes samples from the
    latter halves of multiple independent MCMC chains, then for each chain
    erects a credibility interval and computes the proportion of the
    concatenation of samples from all chains that are covered by the
    credibility interval.  The chains are considered to have converged when the
    mean coverage is within some epsilon of the credibility interval.  This
    convergence condition indicates that the chains are sampling from
    approximately the same distribution.

    References:

      Altekar, G., S. Dwarkadas, J.P. Huelsenbeck, F. Ronquist (2004) Parallel
      Metropolis coupled Markov chain Monte Carlo for Bayesian phylogentic
      inference.  Bioinformatics 20(3):407-415.

      Brooks, S.P., A. Gelman (1998) General methods for monitoring convergence
      of iterative simulations.  J. Comput. Graph. Stat.  7(4):434-455.

      Huelsenbeck, J.P., B. Larget, M.E. Alfaro (2004) Bayesian phylogenetic
      model selection using reversible jump Markov chain Monte Carlo.  Mol.
      Biol. Evol.  21(6):1123-1133.

      Lakner, C., P.v.d. Mark, J.P. Huelsenbeck, B. Larget, F. Ronquist (2008)
      Efficiency of Markov chain Monte Carlo tree proposals in Bayesian
      phylogenetics.  Syst. Biol. 57(1):86-103.

      Lewis, P.O., M.T. Holder, K.E. Holsinger (2005) Polytomies and Bayesian
      phylogenetic inference.  Sys. Biol. 54(2):241-253.
"""

import cPickle
import gc
import os
import random
import sys
import time

import Crux.Config

cdef extern from "Python.h":
    cdef object PyString_FromStringAndSize(char *s, Py_ssize_t len)
from libc cimport *
from libm cimport *
from SFMT cimport *
from Crux.Mc3.Chain cimport *
from Crux.Mc3.Post cimport *
from Crux.Tree cimport Tree, Edge
from Crux.Tree.Bipart cimport Vec, Bipart
from Crux.Tree.Lik cimport Lik
from Crux.Tree.Sumt cimport Part
from Crux.Character cimport Dna

IF @enable_mpi@:
    from mpi4py import MPI
    cimport mpi4py.mpi_c as mpi

    # Use distinct tags for various MPI messaging operations, in order to avoid
    # mixing them up due to out-of-order receipt.
    cdef enum:
        TagHeatSwap = 1
        TagLikLnL   = 2

cdef int _lnLCmp(void *va, void *vb):
    cdef double a = (<double *>va)[0]
    cdef double b = (<double *>vb)[0]

    return (a > b) - (a < b)

# Lookup table of DNA rclasses, used to randomly draw rclasses from the
# appropriate resolution class when drawing a Q from the prior.
cdef list _dnaRclasses = None
cdef void _initRclasses() except *:
    global _dnaRclasses

    assert _dnaRclasses is None
    _dnaRclasses = [
        # 1
        [
            [0,0,0,0,0,0]
        ],
        # 2
        [
            [0,0,0,0,0,1], [0,0,0,0,1,0], [0,0,0,0,1,1], [0,0,0,1,0,0],
            [0,0,0,1,0,1], [0,0,0,1,1,0], [0,0,0,1,1,1], [0,0,1,0,0,0],
            [0,0,1,0,0,1], [0,0,1,0,1,0], [0,0,1,0,1,1], [0,0,1,1,0,0],
            [0,0,1,1,0,1], [0,0,1,1,1,0], [0,0,1,1,1,1], [0,1,0,0,0,0],
            [0,1,0,0,0,1], [0,1,0,0,1,0], [0,1,0,0,1,1], [0,1,0,1,0,0],
            [0,1,0,1,0,1], [0,1,0,1,1,0], [0,1,0,1,1,1], [0,1,1,0,0,0],
            [0,1,1,0,0,1], [0,1,1,0,1,0], [0,1,1,0,1,1], [0,1,1,1,0,0],
            [0,1,1,1,0,1], [0,1,1,1,1,0], [0,1,1,1,1,1]
        ],
        # 3
        [
            [0,0,0,0,1,2], [0,0,0,1,0,2], [0,0,0,1,1,2], [0,0,0,1,2,0],
            [0,0,0,1,2,1], [0,0,0,1,2,2], [0,0,1,0,0,2], [0,0,1,0,1,2],
            [0,0,1,0,2,0], [0,0,1,0,2,1], [0,0,1,0,2,2], [0,0,1,1,0,2],
            [0,0,1,1,1,2], [0,0,1,1,2,0], [0,0,1,1,2,1], [0,0,1,1,2,2],
            [0,0,1,2,0,0], [0,0,1,2,0,1], [0,0,1,2,0,2], [0,0,1,2,1,0],
            [0,0,1,2,1,1], [0,0,1,2,1,2], [0,0,1,2,2,0], [0,0,1,2,2,1],
            [0,0,1,2,2,2], [0,1,0,0,0,2], [0,1,0,0,1,2], [0,1,0,0,2,0],
            [0,1,0,0,2,1], [0,1,0,0,2,2], [0,1,0,1,0,2], [0,1,0,1,1,2],
            [0,1,0,1,2,0], [0,1,0,1,2,1], [0,1,0,1,2,2], [0,1,0,2,0,0],
            [0,1,0,2,0,1], [0,1,0,2,0,2], [0,1,0,2,1,0], [0,1,0,2,1,1],
            [0,1,0,2,1,2], [0,1,0,2,2,0], [0,1,0,2,2,1], [0,1,0,2,2,2],
            [0,1,1,0,0,2], [0,1,1,0,1,2], [0,1,1,0,2,0], [0,1,1,0,2,1],
            [0,1,1,0,2,2], [0,1,1,1,0,2], [0,1,1,1,1,2], [0,1,1,1,2,0],
            [0,1,1,1,2,1], [0,1,1,1,2,2], [0,1,1,2,0,0], [0,1,1,2,0,1],
            [0,1,1,2,0,2], [0,1,1,2,1,0], [0,1,1,2,1,1], [0,1,1,2,1,2],
            [0,1,1,2,2,0], [0,1,1,2,2,1], [0,1,1,2,2,2], [0,1,2,0,0,0],
            [0,1,2,0,0,1], [0,1,2,0,0,2], [0,1,2,0,1,0], [0,1,2,0,1,1],
            [0,1,2,0,1,2], [0,1,2,0,2,0], [0,1,2,0,2,1], [0,1,2,0,2,2],
            [0,1,2,1,0,0], [0,1,2,1,0,1], [0,1,2,1,0,2], [0,1,2,1,1,0],
            [0,1,2,1,1,1], [0,1,2,1,1,2], [0,1,2,1,2,0], [0,1,2,1,2,1],
            [0,1,2,1,2,2], [0,1,2,2,0,0], [0,1,2,2,0,1], [0,1,2,2,0,2],
            [0,1,2,2,1,0], [0,1,2,2,1,1], [0,1,2,2,1,2], [0,1,2,2,2,0],
            [0,1,2,2,2,1], [0,1,2,2,2,2]
        ],
        # 4
        [
            [0,0,0,1,2,3], [0,0,1,0,2,3], [0,0,1,1,2,3], [0,0,1,2,0,3],
            [0,0,1,2,1,3], [0,0,1,2,2,3], [0,0,1,2,3,0], [0,0,1,2,3,1],
            [0,0,1,2,3,2], [0,0,1,2,3,3], [0,1,0,0,2,3], [0,1,0,1,2,3],
            [0,1,0,2,0,3], [0,1,0,2,1,3], [0,1,0,2,2,3], [0,1,0,2,3,0],
            [0,1,0,2,3,1], [0,1,0,2,3,2], [0,1,0,2,3,3], [0,1,1,0,2,3],
            [0,1,1,1,2,3], [0,1,1,2,0,3], [0,1,1,2,1,3], [0,1,1,2,2,3],
            [0,1,1,2,3,0], [0,1,1,2,3,1], [0,1,1,2,3,2], [0,1,1,2,3,3],
            [0,1,2,0,0,3], [0,1,2,0,1,3], [0,1,2,0,2,3], [0,1,2,0,3,0],
            [0,1,2,0,3,1], [0,1,2,0,3,2], [0,1,2,0,3,3], [0,1,2,1,0,3],
            [0,1,2,1,1,3], [0,1,2,1,2,3], [0,1,2,1,3,0], [0,1,2,1,3,1],
            [0,1,2,1,3,2], [0,1,2,1,3,3], [0,1,2,2,0,3], [0,1,2,2,1,3],
            [0,1,2,2,2,3], [0,1,2,2,3,0], [0,1,2,2,3,1], [0,1,2,2,3,2],
            [0,1,2,2,3,3], [0,1,2,3,0,0], [0,1,2,3,0,1], [0,1,2,3,0,2],
            [0,1,2,3,0,3], [0,1,2,3,1,0], [0,1,2,3,1,1], [0,1,2,3,1,2],
            [0,1,2,3,1,3], [0,1,2,3,2,0], [0,1,2,3,2,1], [0,1,2,3,2,2],
            [0,1,2,3,2,3], [0,1,2,3,3,0], [0,1,2,3,3,1], [0,1,2,3,3,2],
            [0,1,2,3,3,3]
        ],
        # 5
        [
            [0,0,1,2,3,4], [0,1,0,2,3,4], [0,1,1,2,3,4], [0,1,2,0,3,4],
            [0,1,2,1,3,4], [0,1,2,2,3,4], [0,1,2,3,0,4], [0,1,2,3,1,4],
            [0,1,2,3,2,4], [0,1,2,3,3,4], [0,1,2,3,4,0], [0,1,2,3,4,1],
            [0,1,2,3,4,2], [0,1,2,3,4,3], [0,1,2,3,4,4]
        ],
        # 6
        [
            [0,1,2,3,4,5]
        ]
    ]

cdef class Mc3:
    """
        alignment
          Aligned character-by-taxon matrix input data.

        outPrefix
          Filename prefix, used for naming output files.
    """
    def __cinit__(self):
        cdef unsigned i

        self.swapInfo = NULL
        self.swapStats = NULL
        for 0 <= i < PropCnt:
            self.propStats[i] = NULL
        self.cachedLnLs = NULL
        self.lnLs = NULL
        IF @enable_mpi@:
            self.mpiLeaderCommAlloced = False
            self.mpiChainComms = NULL

    def __dealloc__(self):
        cdef unsigned i

        if self.swapInfo != NULL:
            free(self.swapInfo)
            self.swapInfo = NULL
        if self.swapStats != NULL:
            free(self.swapStats)
            self.swapStats = NULL
        for 0 <= i < PropCnt:
            if self.propStats[i] != NULL:
                free(self.propStats[i])
                self.propStats[i] = NULL
        if self.cachedLnLs != NULL:
            free(self.cachedLnLs)
            self.cachedLnLs = NULL
        if self.lnLs != NULL:
            for 0 <= i < self._nruns+1:
                if self.lnLs[i] != NULL:
                    free(self.lnLs[i])
                    self.lnLs[i] = NULL
            free(self.lnLs)
            self.lnLs = NULL
        IF @enable_mpi@:
            if self.mpiLeaderCommAlloced:
                mpi.MPI_Comm_free(&self.mpiLeaderComm)
                self.mpiLeaderCommAlloced = False
            if self.mpiChainComms != NULL:
                for 0 <= i < self._nruns * self._ncoupled:
                    mpi.MPI_Comm_free(&self.mpiChainComms[i])
                free(self.mpiChainComms)
                self.mpiChainComms = NULL

    def __init__(self, Alignment alignment, str outPrefix):
        self.alignment = alignment
        self.outPrefix = outPrefix

        # Set defaults.
        self._graphDelay = -1.0
        self._emaAlpha = 1.0
        self._cvgSampStride = 1
        self._cvgAlpha = 0.05
        self._cvgEpsilon = 0.01
        self._minStep = 100000
        self._maxStep = ULLONG_MAX
        self._stride = 100
        self._nruns = 2
        self._ncoupled = 1
        self._heatDelta = 0.05
        self._swapStride = 1
        self._nmodels = 1
        self._ncat = 1
        self._catMedian = False
        self._invar = False
        self._weightLambda = 2.0 * log(1.6)
        self._freqLambda = 2.0 * log(1.6)
        self._rmultLambda = 2.0 * log(1.6)
        self._rateLambda = 2.0 * log(1.6)
        self._rateShapeInvLambda = 2.0 * log(1.6)
        self._invarLambda = 2.0 * log(1.6)
        self._brlenLambda = 2.0 * log(1.6)
        self._etbrPExt = 0.8
        self._etbrLambda = 2.0 * log(1.6)
        self._rateShapeInvPrior = 1.0
        self._invarPrior = 0.5
        self._brlenPrior = 10.0
        self._rateJumpPrior = 1.0
        self._polytomyJumpPrior = 1.0
        self._rateShapeInvJumpPrior = 1.0
        self._invarJumpPrior = 1.0
        self._freqJumpPrior = 1.0
        self._mixtureJumpPrior = 1.0 / 3.0
        self.props[PropWeight] = 30.0
        self.props[PropFreq] = 10.0
        self.props[PropRmult] = 10.0
        self.props[PropRate] = 30.0
        self.props[PropRateShapeInv] = 10.0
        self.props[PropInvar] = 5.0
        self.props[PropBrlen] = 20.0
        self.props[PropEtbr] = 80.0
        self.props[PropRateJump] = 10.0
        self.props[PropPolytomyJump] = 40.0
        self.props[PropRateShapeInvJump] = 10.0
        self.props[PropInvarJump] = 5.0
        self.props[PropFreqJump] = 10.0
        self.props[PropMixtureJump] = 5.0

    cpdef Mc3 dup(self):
        """
            Create a duplicate instance that encapsulates parameters, but
            excludes the posterior distribution.
        """
        cdef Mc3 ret
        cdef unsigned i

        ret = Mc3(self.alignment, self.outPrefix)
        ret._graphDelay = self._graphDelay
        ret._emaAlpha = self._emaAlpha
        ret._cvgSampStride = self._cvgSampStride
        ret._cvgAlpha = self._cvgAlpha
        ret._cvgEpsilon = self._cvgEpsilon
        ret._minStep = self._minStep
        ret._maxStep = self._maxStep
        ret._stride = self._stride
        ret._nruns = self._nruns
        ret._ncoupled = self._ncoupled
        ret._heatDelta = self._heatDelta
        ret._swapStride = self._swapStride
        ret._nmodels = self._nmodels
        ret._ncat = self._ncat
        ret._catMedian = self._catMedian
        ret._invar = self._invar
        ret._weightLambda = self._weightLambda
        ret._freqLambda = self._freqLambda
        ret._rmultLambda = self._rmultLambda
        ret._rateLambda = self._rateLambda
        ret._rateShapeInvLambda = self._rateShapeInvLambda
        ret._invarLambda = self._invarLambda
        ret._brlenLambda = self._brlenLambda
        ret._etbrPExt = self._etbrPExt
        ret._etbrLambda = self._etbrLambda
        ret._rateShapeInvPrior = self._rateShapeInvPrior
        ret._invarPrior = self._invarPrior
        ret._brlenPrior = self._brlenPrior
        ret._rateJumpPrior = self._rateJumpPrior
        ret._polytomyJumpPrior = self._polytomyJumpPrior
        ret._rateShapeInvJumpPrior = self._rateShapeInvJumpPrior
        ret._invarJumpPrior = self._invarJumpPrior
        ret._freqJumpPrior = self._freqJumpPrior
        ret._mixtureJumpPrior = self._mixtureJumpPrior

        for 0 <= i < PropCnt:
            ret.props[i] = self.props[i]

    IF @enable_mpi@:
        cdef void storeSwapStats(self) except *:
            cdef unsigned i
            cdef Mc3RateStats *swapStats
            cdef uint64_t n

            if self.mpiLeaderRank >= 0:
                for 0 <= i < self._nruns:
                    swapStats = &self.swapStats[i]
                    n = swapStats.n
                    mpi.MPI_Reduce(&n, &swapStats.n, 1, \
                      mpi.MPI_UNSIGNED_LONG_LONG, mpi.MPI_SUM, 0, \
                      self.mpiLeaderComm)

        cdef void storePropStats(self) except *:
            cdef unsigned i, j
            cdef Mc3RateStats *propStats
            cdef uint64_t nd[2]

            if self.mpiLeaderRank >= 0:
                for 0 <= i < PropCnt:
                    for 0 <= j < self._nruns:
                        propStats = &self.propStats[i][j]
                        nd[0] = propStats.n
                        nd[1] = propStats.d
                        mpi.MPI_Reduce(nd, &propStats.n, 2, \
                          mpi.MPI_UNSIGNED_LONG_LONG, mpi.MPI_SUM, 0, \
                          self.mpiLeaderComm)

    # Update diagnostics based on swapStats and propStats.
    cdef void updateDiags(self, uint64_t step) except *:
        cdef unsigned runInd, i
        cdef Mc3RateStats *rateStats
        cdef double xt

        IF @enable_mpi@:
            if self.mpiLeaderRank != 0:
                # Only the master node has all the information to compute the
                # diagnostics.  This node only needs to clear its partial
                # stats.
                for 0 <= runInd < self._nruns:
                    rateStats = &self.swapStats[runInd]
                    rateStats.n = 0
                    for 0 <= i < PropCnt:
                        rateStats = &self.propStats[i][runInd]
                        rateStats.n = 0
                        rateStats.d = 0
                return

        for 0 <= runInd < self._nruns:
            # Compute the exponential moving averages for heat swap rates.
            # Swaps are double-counted, since both participants count them.
            rateStats = &self.swapStats[runInd]
            rateStats.ema += self._emaAlpha * \
              (<double>rateStats.n/(2.0*<double>self._stride) - \
              rateStats.ema)
            rateStats.n = 0

            # Compute the exponential moving averages for proposal acceptance
            # rates.
            for 0 <= i < PropCnt:
                rateStats = &self.propStats[i][runInd]
                if rateStats.d == 0:
                    xt = 0.0
                else:
                    xt = (<double>rateStats.n/<double>rateStats.d)
                rateStats.ema += self._emaAlpha * (xt - rateStats.ema)
                IF @enable_mpi@:
                    rateStats.n = 0
                    rateStats.d = 0

    cdef void storeLiksLnLsUni(self, uint64_t step) except *:
        cdef uint64_t samp
        cdef unsigned i, j
        cdef size_t lnLsMax
        cdef double *lnLs

        samp = step / self._stride

        # Allocate/expand lnLs as necessary.
        if self.lnLsMax == 0:
            if self._minStep == 0:
                lnLsMax = 8 # Arbitrary starting size.
            else:
                lnLsMax = self._minStep / self._stride
            for 0 <= i < self._nruns+1:
                self.lnLs[i] = <double *>malloc(lnLsMax * sizeof(double))
                if self.lnLs[i] == NULL:
                    for 0 <= j < i:
                        free(self.lnLs[j])
                    raise MemoryError("Error allocating lnLs[%d]" % i)
            self.lnLsMax = lnLsMax
        elif self.lnLsMax <= samp:
            lnLsMax = self.lnLsMax + 4 + (self.lnLsMax >> 3) + \
              (self.lnLsMax >> 4)
            for 0 <= i < self._nruns+1:
                lnLs = <double *>realloc(self.lnLs[i], lnLsMax * sizeof(double))
                if lnLs == NULL:
                    raise MemoryError("Error reallocating lnLs[%d]" % i)
                self.lnLs[i] = lnLs
            self.lnLsMax = lnLsMax

        # Store lnLs.
        for 0 <= i < self._nruns:
            self.lnLs[i][samp] = self.cachedLnLs[i]

    IF @enable_mpi@:
        cdef void storeLiksLnLsMpi(self, uint64_t step) except *:
            cdef unsigned i, nRecv, iRecv
            cdef int pSize
            cdef Lik lik, masterLik
            cdef mpi.MPI_Status status
            cdef str pickle
            cdef char *s
            cdef object tuple
            IF @enable_debug@:
                cdef double lnL

            if self.mpiLeaderRank == -1:
                return

            # Force synchronization between the master and slaves.  Without
            # this synchronization, it would be possible for messages from
            # different steps to be interleaved.
            mpi.MPI_Barrier(self.mpiLeaderComm)

            if self.mpiLeaderRank != 0:
                for 0 <= i < self._nruns:
                    # If a node other than the master owns the unheated chain
                    # for run i, transmit the lik/lnL to the master.  Since
                    # heats are swapped among coupled chains, there's no way to
                    # know a priori which chain is unheated, so start off by
                    # determining the owner of the unheated chain.
                    lik = <Lik>self.liks[i]
                    if lik is not None:
                        # Send lik/lnL to the master.
                        tuple = (i, cPickle.dumps(lik), self.cachedLnLs[i])
                        IF @enable_debug@:
                            lnL = lik.unpickle(tuple[1]).lnL()
                            assert (0.99 < self.cachedLnLs[i]/lnL) and \
                              (self.cachedLnLs[i]/lnL < 1.01)
                        pickle = cPickle.dumps(tuple)
                        s = pickle
                        pSize = len(pickle)
                        mpi.MPI_Send(s, pSize, mpi.MPI_BYTE, 0, TagLikLnL, \
                          self.mpiLeaderComm)
            else:
                # Obtain a reference to a Lik that was constructed using the
                # same alignment as those that will be unpickled below.
                masterLik = (<Chain>(<list>self.runs[0])[0]).lik
                assert masterLik is not None

                # Count how many messages to receive.
                nRecv = 0
                for 0 <= i < self._nruns:
                    if self.liks[i] is None:
                        nRecv += 1

                for 0 <= i < nRecv:
                    # Receive lik/lnL from a slave.  Messages may arrive out of
                    # order, so use iRecv rather linear indexing.
                    assert self.mpiLeaderRank == 0

                    # Get pickle size.
                    mpi.MPI_Probe(mpi.MPI_ANY_SOURCE, TagLikLnL, \
                      self.mpiLeaderComm, &status)
                    mpi.MPI_Get_count(&status, mpi.MPI_BYTE, &pSize)
                    s = <char *>malloc(pSize)
                    if s == NULL:
                        raise MemoryError("Error allocating Lik pickle string")
                    # Get pickle.
                    mpi.MPI_Recv(s, pSize, mpi.MPI_BYTE, status.MPI_SOURCE, \
                      TagLikLnL, self.mpiLeaderComm, mpi.MPI_STATUS_IGNORE)
                    pickle = PyString_FromStringAndSize(s, pSize)
                    free(s)
                    # Unpickle.
                    tuple = cPickle.loads(pickle)
                    iRecv = <unsigned>tuple[0]
                    assert self.liks[iRecv] is None
                    self.liks[iRecv] = masterLik.unpickle(<str>tuple[1])
                    self.liks[iRecv].prep()
                    self.cachedLnLs[iRecv] = <double>tuple[2]
                    IF @enable_debug@:
                        lnL = self.liks[iRecv].lnL()
                        assert (0.99 < self.cachedLnLs[iRecv]/lnL) and \
                          (self.cachedLnLs[iRecv]/lnL < 1.01)

                IF @enable_debug@:
                    for 0 <= i < self._nruns:
                        assert self.liks[i] is not None
                self.storeLiksLnLsUni(step)

    cdef void storeLiksLnLs(self, uint64_t step) except *:
        IF @enable_mpi@:
            if self.mpiWorldSize == 1:
                self.storeLiksLnLsUni(step)
            else:
                self.storeLiksLnLsMpi(step)
        ELSE:
            self.storeLiksLnLsUni(step)

    cdef void sendSample(self, unsigned runInd, uint64_t step, \
      double heat, uint64_t nswap, uint64_t *accepts, uint64_t *rejects, \
      Lik lik, double lnL) except *:
        cdef unsigned i

        # Record swap statistics.
        self.swapStats[runInd].n += nswap

        if heat != 1.0:
            # The remainder of this method applies only to cold chains.
            return

        # Record accept/reject rate statistics.
        for 0 <= i < PropCnt:
            self.propStats[i][runInd].n = accepts[i]
            self.propStats[i][runInd].d = accepts[i] + rejects[i]

        self.liks[runInd] = lik

        self.cachedLnLs[runInd] = lnL

    cdef void sendSwapInfo(self, unsigned runInd, unsigned srcChainInd, \
      unsigned dstChainInd, uint64_t step, double heat, double lnL) except *:
        cdef Mc3SwapInfo *swapInfo

        assert step % self._swapStride == 0
        swapInfo = &self.swapInfo[((runInd*self._ncoupled + \
          srcChainInd)*self._ncoupled + dstChainInd)*2 + \
          ((step/self._swapStride) % 2)]
        assert swapInfo.step == 0

        swapInfo.step = step
        swapInfo.heat = heat
        swapInfo.lnL = lnL

    cdef void recvSwapInfoUni(self, unsigned runInd, unsigned dstChainInd, \
      unsigned srcChainInd, uint64_t step, double *heat, double *lnL) except *:
        cdef Mc3SwapInfo *swapInfo

        assert step % self._swapStride == 0
        swapInfo = &self.swapInfo[((runInd*self._ncoupled + \
          srcChainInd)*self._ncoupled + dstChainInd)*2 + \
          ((step/self._swapStride) % 2)]
        assert swapInfo.step == step

        heat[0] = swapInfo.heat
        lnL[0] = swapInfo.lnL

        IF @enable_debug@:
            swapInfo.step = 0

    IF @enable_mpi@:
        cdef void recvSwapInfoMpi(self, unsigned runInd, unsigned dstChainInd, \
          unsigned srcChainInd, uint64_t step, double *heat, double *lnL) \
          except *:
            cdef Mc3SwapInfo *swapInfoSend, *swapInfoRecv
            cdef int peerRank
            cdef unsigned chain

            swapInfoSend = &self.swapInfo[((runInd*self._ncoupled + \
              dstChainInd)*self._ncoupled + srcChainInd)*2 + \
              ((step/self._swapStride) % 2)]

            swapInfoRecv = &self.swapInfo[((runInd*self._ncoupled + \
              srcChainInd)*self._ncoupled + dstChainInd)*2 + \
              ((step/self._swapStride) % 2)]

            peerRank = (runInd*self._ncoupled + srcChainInd) % \
              self.mpiLeaderSize
            if peerRank != (runInd*self._ncoupled + dstChainInd) % \
              self.mpiLeaderSize:
                assert swapInfoSend.step == step
                if self.mpiLeaderRank >= 0:
                    assert (runInd*self._ncoupled + dstChainInd) % \
                      self.mpiLeaderSize == self.mpiLeaderRank
                    mpi.MPI_Sendrecv(swapInfoSend, sizeof(Mc3SwapInfo), \
                      mpi.MPI_BYTE, peerRank, TagHeatSwap, swapInfoRecv, \
                      sizeof(Mc3SwapInfo), mpi.MPI_BYTE, peerRank, \
                      TagHeatSwap, self.mpiLeaderComm, mpi.MPI_STATUS_IGNORE)
                IF @enable_debug@:
                    swapInfoSend.step = 0

            # Share swap info with follower nodes.
            chain = runInd*self._ncoupled + dstChainInd
            mpi.MPI_Bcast(swapInfoRecv, sizeof(Mc3SwapInfo), mpi.MPI_BYTE,
              0, self.mpiChainComms[chain])

            self.recvSwapInfoUni(runInd, dstChainInd, srcChainInd, step, \
              heat, lnL)

    cdef void recvSwapInfo(self, unsigned runInd, unsigned dstChainInd, \
      unsigned srcChainInd, uint64_t step, double *heat, double *lnL) except *:
        IF @enable_mpi@:
            if self.mpiWorldSize == 1:
                self.recvSwapInfoUni(runInd, dstChainInd, srcChainInd, step, \
                  heat, lnL)
            else:
                self.recvSwapInfoMpi(runInd, dstChainInd, srcChainInd, step, \
                  heat, lnL)
        ELSE:
            self.recvSwapInfoUni(runInd, dstChainInd, srcChainInd, step, heat, \
              lnL)

    IF @enable_mpi@:
        cdef void initComms(self) except *:
            cdef unsigned r, i, j, k, nchains, chain
            cdef list leaderList, chainLists
            cdef int *leaderArray, *activeArray, *chainArray
            cdef mpi.MPI_Group worldGroup, leaderGroup, activeGroup, chainGroup

            # Store the number of nodes, and this node's rank.
            mpi.MPI_Comm_size(mpi.MPI_COMM_WORLD, &self.mpiWorldSize)
            mpi.MPI_Comm_rank(mpi.MPI_COMM_WORLD, &self.mpiWorldRank)

            # Create a list of chain leaders that is ordered such that chains
            # will be interleaved if there are more chains than nodes.
            leaderList = []
            chainLists = []
            for 0 <= i < self._nruns:
                for 0 <= j < self._ncoupled:
                    r = (i*self._ncoupled + j) % self.mpiWorldSize
                    if r not in leaderList:
                        leaderList.append(r)
                    chainLists.append([])

            # Create lists of nodes associated with chains.
            nchains = self._nruns * self._ncoupled
            if self.mpiWorldSize <= nchains:
                r = 0
                for 0 <= i < self._nruns:
                    for 0 <= j < self._ncoupled:
                        chain = i*self._ncoupled + j
                        assert r not in chainLists[chain]
                        chainLists[chain].append(r)
                        r = (r + 1) % self.mpiWorldSize
            else:
                for 0 <= r < self.mpiWorldSize:
                    chain = r % nchains
                    assert r not in chainLists[chain]
                    chainLists[chain].append(r)

            mpi.MPI_Comm_group(mpi.MPI_COMM_WORLD, &worldGroup)
            # Create the leader communicator.
            leaderArray = <int *>malloc(len(leaderList) * sizeof(int))
            if leaderArray == NULL:
                raise MemoryError("Error allocating leaderArray")
            for 0 <= i < len(leaderList):
                leaderArray[i] = <int>leaderList[i]

            try:
                mpi.MPI_Group_incl(worldGroup, len(leaderList), \
                  leaderArray, &leaderGroup)
                mpi.MPI_Comm_create(mpi.MPI_COMM_WORLD, leaderGroup, \
                  &self.mpiLeaderComm)
                self.mpiLeaderCommAlloced = True
                mpi.MPI_Group_free(&leaderGroup)
            finally:
                free(leaderArray)

            # Store the number of leader nodes, and this node's rank.
            self.mpiLeaderSize = len(leaderList)
            if self.mpiWorldRank in leaderList:
                mpi.MPI_Comm_rank(self.mpiLeaderComm, &self.mpiLeaderRank)
            else:
                self.mpiLeaderRank = -1

            # Create a communicator for each chain.
            self.mpiChainComms = <mpi.MPI_Comm *>malloc(nchains * \
              sizeof(mpi.MPI_Comm))
            if self.mpiChainComms == NULL:
                raise MemoryError("Error allocating mpiChainComms")

            try:
                for 0 <= i < self._nruns:
                    for 0 <= j < self._ncoupled:
                        chain = i*self._ncoupled + j
                        chainArray = <int *>malloc(len(chainLists[chain]) * \
                          sizeof(int))
                        if chainArray == NULL:
                            raise MemoryError( \
                              "Error allocating chainArray")
                        for 0 <= k < len(chainLists[chain]):
                            chainArray[k] = <int>chainLists[chain][k]

                        try:
                            mpi.MPI_Group_incl(worldGroup, \
                              len(chainLists[chain]), chainArray, &chainGroup)
                            mpi.MPI_Comm_create(mpi.MPI_COMM_WORLD, \
                              chainGroup, &self.mpiChainComms[chain])
                            mpi.MPI_Group_free(&chainGroup)
                        finally:
                            free(chainArray)
            except:
                for 0 <= i < chain:
                    mpi.MPI_Comm_free(&self.mpiChainComms[i])
                free(self.mpiChainComms)
                self.mpiChainComms = NULL
                raise

    cdef void initLogs(self) except *:
        cdef file f
        cdef unsigned nchars, j

        IF @enable_mpi@:
            if self.mpiLeaderRank != 0:
                return

        nchars = 0
        for 0 <= j < self.alignment.nchars:
            nchars += self.alignment.getFreq(j)

        self.lFile = open("%s.l" % self.outPrefix, "w")
        for f in ((self.lFile, sys.stdout) if self.verbose else (self.lFile,)):
            f.write("Begin run: %s\n" % \
              time.strftime("%Y/%m/%d %H:%M:%S (%Z)", \
              time.localtime(time.time())))
            f.write("Crux version: @crux_version@\n")
            f.write("Host machine: %r\n" % (os.uname(),))
            IF @enable_mpi@:
                f.write("MPI node count: %d\n" % self.mpiWorldSize)
            f.write("PRNG seed: %d\n" % Crux.Config.seed)
            f.write("Taxa: %d\n" % self.alignment.ntaxa)
            f.write("Characters: %d\n" % nchars)
            f.write("Unique site patterns: %d\n" % self.alignment.nchars)

            f.write("Configuration parameters:\n")
            f.write("  outPrefix: %r\n" % self.outPrefix)
            f.write("  graphDelay: %r\n" % self._graphDelay)
            f.write("  emaAlpha: %r\n" % self._emaAlpha)
            f.write("  cvgSampStride: %r\n" % self._cvgSampStride)
            f.write("  cvgAlpha: %r\n" % self._cvgAlpha)
            f.write("  cvgEpsilon: %r\n" % self._cvgEpsilon)
            f.write("  minStep: %r\n" % self._minStep)
            f.write("  maxStep: %r\n" % self._maxStep)
            f.write("  stride: %r\n" % self._stride)
            f.write("  nruns: %r\n" % self._nruns)
            f.write("  ncoupled: %r\n" % self._ncoupled)
            f.write("  heatDelta: %r\n" % self._heatDelta)
            f.write("  swapStride: %r\n" % self._swapStride)
            f.write("  nmodels: %r\n" % self._nmodels)
            f.write("  ncat: %r\n" % self._ncat)
            f.write("  catMedian: %r\n" % self._catMedian)
            f.write("  invar: %r\n" % self._invar)
            f.write("  weightLambda: %r\n" % self._weightLambda)
            f.write("  freqLambda: %r\n" % self._freqLambda)
            f.write("  rmultLambda: %r\n" % self._rmultLambda)
            f.write("  rateLambda: %r\n" % self._rateLambda)
            f.write("  rateShapeInvLambda: %r\n" % self._rateShapeInvLambda)
            f.write("  invarLambda: %r\n" % self._invarLambda)
            f.write("  brlenLambda: %r\n" % self._brlenLambda)
            f.write("  etbrPExt: %r\n" % self._etbrPExt)
            f.write("  etbrLambda: %r\n" % self._etbrLambda)
            f.write("  rateShapeInvPrior: %r\n" % self._rateShapeInvPrior)
            f.write("  invarPrior: %r\n" % self._invarPrior)
            f.write("  brlenPrior: %r\n" % self._brlenPrior)
            f.write("  rateJumpPrior: %r\n" % self._rateJumpPrior)
            f.write("  polytomyJumpPrior: %r\n" % self._polytomyJumpPrior)
            f.write("  rateShapeInvJumpPrior: %r\n" % \
              self._rateShapeInvJumpPrior)
            f.write("  invarJumpPrior: %r\n" % self._invarJumpPrior)
            f.write("  freqJumpPrior: %r\n" % self._freqJumpPrior)
            f.write("  mixtureJumpPrior: %r\n" % self._mixtureJumpPrior)
            f.write("  weightProp: %r\n" % self.props[PropWeight])
            f.write("  freqProp: %r\n" % self.props[PropFreq])
            f.write("  rmultProp: %r\n" % self.props[PropRmult])
            f.write("  rateProp: %r\n" % self.props[PropRate])
            f.write("  rateShapeInvProp: %r\n" % \
              self.props[PropRateShapeInv])
            f.write("  invarProp: %r\n" % self.props[PropInvar])
            f.write("  brlenProp: %r\n" % self.props[PropBrlen])
            f.write("  etbrProp: %r\n" % self.props[PropEtbr])
            f.write("  rateJumpProp: %r\n" % self.props[PropRateJump])
            f.write("  polytomyJumpProp: %r\n" % self.props[PropPolytomyJump])
            f.write("  rateShapeInvJumpProp: %r\n" % \
              self.props[PropRateShapeInvJump])
            f.write("  invarJumpProp: %r\n" % self.props[PropInvarJump])
            f.write("  freqJumpProp: %r\n" % self.props[PropFreqJump])
            f.write("  mixtureJumpProp: %r\n" % self.props[PropMixtureJump])
        self.lFile.flush()

        self.tFile = open("%s.t" % self.outPrefix, "w")
        self.tFile.write("[[run step] tree]\n")
        self.tFile.flush()

        self.pFile = open("%s.p" % self.outPrefix, "w")
        self.pFile.write( \
          "run\tstep\tmodel\tweight\trmult\t( rclass )[ R ] wNorm alpha" \
          " pinvar [ Pi ]\n")
        self.pFile.flush()
        if self.verbose:
            sys.stdout.write( \
              "p\trun\tstep\tmodel\tweight\trmult\t( rclass )[ R ] wNorm" \
              " alpha pinvar [ Pi ]\n")

        self.sFile = open("%s.s" % self.outPrefix, "w")
        self.sFile.write("step\t[ lnLs ] Rcov [ swapRates ] {" \
          " [ weightPropRates ] [ freqPropRates ] [ rmultPropRates ]" \
          " [ ratePropRates ] [ rateShapeInvPropRates ] [ invarPropRates ]" \
          " [ brlenPropRates ] [ etbrPropRates ] [ rateJumpPropRates ]" \
          " [ polytomyJumpPropRates ] [ rateShapeInvJumpPropRates ]" \
          " [ invarJumpPropRates ] [ freqJumpPropRates ]" \
          " [ mixtureJumpPropRates ] }\n")
        self.sFile.flush()
        if self.verbose:
            sys.stdout.write("s\tstep\t[ lnLs ] Rcov\n")

    # Allocate swapInfo matrix.
    cdef void initSwapInfo(self) except *:
        if self._ncoupled > 1:
            if self.swapInfo != NULL:
                free(self.swapInfo)
            self.swapInfo = <Mc3SwapInfo *>calloc(self._nruns * \
              self._ncoupled * self._ncoupled * 2, sizeof(Mc3SwapInfo))
            if self.swapInfo == NULL:
                raise MemoryError("Error allocating swapInfo")

    # Allocate swapStats vector.
    cdef void initSwapStats(self) except *:
        if self.swapStats != NULL:
            free(self.swapStats)
        self.swapStats = <Mc3RateStats *>calloc(self._nruns, \
          sizeof(Mc3RateStats))
        if self.swapStats == NULL:
            raise MemoryError("Error allocating swapStats")

    # Allocate propStats vectors.
    cdef void initPropStats(self) except *:
        cdef unsigned i

        for 0 <= i < PropCnt:
            if self.propStats[i] != NULL:
                free(self.propStats[i])
            self.propStats[i] = <Mc3RateStats *>calloc(self._nruns, \
              sizeof(Mc3RateStats))
            if self.propStats[i] == NULL:
                raise MemoryError("Error allocating propStats[%d]" % i)

    cdef void initLiks(self) except *:
        self.liks = [None] * self._nruns

    cdef void resetLiks(self) except *:
        cdef unsigned i

        for 0 <= i < self._nruns:
            self.liks[i] = None

    cdef void initLnLs(self) except *:
        if self.cachedLnLs != NULL:
            free(self.cachedLnLs)
            self.cachedLnLs = NULL
        if self.lnLs != NULL:
            free(self.lnLs)
            self.lnLs = NULL
        self.lnLsMax = 0

        self.cachedLnLs = <double *>malloc(self._nruns * sizeof(double))
        if self.cachedLnLs == NULL:
            raise MemoryError("Error allocating cachedLnLs")

        self.lnLs = <double **>calloc(self._nruns+1, sizeof(double *))
        if self.lnLs == NULL:
            raise MemoryError("Error allocating lnLs")

    # Initialize props[PC]df, which are used to compute proposal ratios and
    # choose proposal types.
    cdef void initProps(self) except *:
        cdef double propsSum
        cdef unsigned i

        assert sizeof(self.props)/sizeof(double) == PropCnt
        assert sizeof(self.propsPdf)/sizeof(double) == PropCnt
        assert sizeof(self.propsCdf)/sizeof(double) == PropCnt + 1
        if self.alignment.charType is not Dna:
            # Mixture jumps are only implemented for DNA.
            self.props[PropMixtureJump] = 0.0
        if self._nmodels < 2 and self.props[PropMixtureJump] == 0.0:
            # There are not enough models in the mixture for weights or rate
            # multipliers to be relevant.
            self.props[PropWeight] = 0.0
            self.props[PropRmult] = 0.0
        if self.alignment.charType.get().nstates <= 2:
            # There are not enough states to allow rate class grouping.
            self.props[PropRateJump] = 0.0
        if self.alignment.ntaxa < 3:
            # There are not enough taxa to allow multiple simultaneous branch
            # length changes.
            self.props[PropEtbr] = 0.0
        if self.alignment.ntaxa < 4:
            # There are not enough taxa to allow topology changes.
            self.props[PropPolytomyJump] = 0.0
        if self._ncat < 2:
            # There aren't multiple rate categories, so Gamma-distributed rates
            # are irrelevant.
            self.props[PropRateShapeInv] = 0.0
            self.props[PropRateShapeInvJump] = 0.0
        if not self._invar:
            # Invariable sites are not supported.
            self.props[PropInvar] = 0.0
            self.props[PropInvarJump] = 0.0
        propsSum = 0.0
        for 0 <= i < PropCnt:
            propsSum += self.props[i]
        if propsSum == 0.0:
            raise ValueError("No proposals are enabled")
        for 0 <= i < PropCnt:
            self.propsPdf[i] = self.props[i] / propsSum
        self.propsCdf[0] = 0.0
        for 1 <= i < PropCnt+1:
            self.propsCdf[i] = self.propsCdf[i-1] + (self.props[i-1] / propsSum)

    cdef void initRunsUni(self, list liks) except *:
        cdef uint32_t seed
        cdef unsigned i, j, chain
        cdef list swapSeeds, chainSeeds, run

        seed = random.randint(0, 0xffffffffU)
        # Create a swap seed for each set of coupled chains.
        swapSeeds = []
        for 0 <= i < self._nruns:
            swapSeeds.append(seed)
            seed += 1
        chainSeeds = []
        # Create one seed for each set of nodes that handle a chain.
        nchains = self._nruns * self._ncoupled
        for 0 <= i < nchains:
            chainSeeds.append(seed)
            seed += 1

        # Create a list of lists that will contain references to all chains.
        self.runs = []
        for 0 <= i < self._nruns:
            run = []
            self.runs.append(run)
            for 0 <= j < self._ncoupled:
                run.append(None)

        # Initialize all chains.
        for 0 <= i < self._nruns:
            for 0 <= j < self._ncoupled:
                chain = i*self._ncoupled + j
                # Create the chain.
                self.runs[i][j] = Chain(self, i, j, swapSeeds[i], \
                  chainSeeds[chain], liks[chain])

    IF @enable_mpi@:
        cdef void initRunsMpi(self, list liks) except *:
            cdef uint32_t seed
            cdef unsigned i, j, nchains, r, chain
            cdef list swapSeeds, chainSeeds, run

            seed = random.randint(0, 0xffffffffU)
            # Create a swap seed for each set of coupled chains.
            swapSeeds = []
            for 0 <= i < self._nruns:
                swapSeeds.append(seed)
                seed += 1
            chainSeeds = []
            # Create one seed for each set of nodes that handle a chain.
            nchains = self._nruns * self._ncoupled
            for 0 <= i < nchains:
                chainSeeds.append(seed)
                seed += 1

            # Create a list of lists that will contain references to any chains
            # this node handles.
            self.runs = []
            for 0 <= i < self._nruns:
                run = []
                self.runs.append(run)
                for 0 <= j < self._ncoupled:
                    run.append(None)

            # Initialize the chains this node handles.
            if self.mpiWorldSize <= nchains:
                for 0 <= i < self._nruns:
                    for 0 <= j < self._ncoupled:
                        chain = i*self._ncoupled + j
                        if chain % self.mpiWorldSize == self.mpiWorldRank:
                            # Configure the Lik for MPI striping.
                            (<Lik>liks[chain]).configMpi( \
                              self.mpiChainComms[chain])
                            # Create the chain.
                            self.runs[i][j] = Chain(self, i, j, swapSeeds[i], \
                              chainSeeds[chain], liks[chain])
            else:
                r = self.mpiWorldRank
                i = (r % nchains) / self._ncoupled
                j = (r % nchains) % self._ncoupled
                chain = i*self._ncoupled + j
                # Configure the Lik for MPI striping.
                (<Lik>liks[r % nchains]).configMpi(self.mpiChainComms[chain])
                # Create the chain.
                self.runs[i][j] = Chain(self, i, j, swapSeeds[i], \
                  chainSeeds[r % nchains], liks[r % nchains])

    # Create/initialize chains.
    cdef void initRuns(self, list liks) except *:
        cdef unsigned i

        if liks is None:
            liks = [self.randomLik() for i in \
              xrange(self._nruns * self._ncoupled)]
        assert len(liks) == self._nruns * self._ncoupled

        IF @enable_mpi@:
            if self.mpiWorldSize == 1:
                self.initRunsUni(liks)
            else:
                self.initRunsMpi(liks)
        ELSE:
            self.initRunsUni(liks)

    cdef void advanceUni(self) except *:
        cdef unsigned i, j
        cdef list run
        cdef Chain chain

        for 0 <= i < self._nruns:
            run = <list>self.runs[i]
            for 0 <= j < self._ncoupled:
                chain = <Chain>run[j]
                chain.advance0()

        for 0 <= i < self._nruns:
            run = <list>self.runs[i]
            for 0 <= j < self._ncoupled:
                chain = <Chain>run[j]
                chain.advance1()

    IF @enable_mpi@:
        cdef void advanceMpi(self) except *:
            cdef unsigned nchains, i, j, r
            cdef list run
            cdef Chain chain

            nchains = self._nruns * self._ncoupled
            if self.mpiWorldSize <= nchains:
                for 0 <= i < self._nruns:
                    for 0 <= j < self._ncoupled:
                        if (i*self._ncoupled + j) % self.mpiWorldSize == \
                          self.mpiWorldRank:
                            run = <list>self.runs[i]
                            chain = <Chain>run[j]
                            chain.advance0()
                for 0 <= i < self._nruns:
                    for 0 <= j < self._ncoupled:
                        if (i*self._ncoupled + j) % self.mpiWorldSize == \
                          self.mpiWorldRank:
                            run = <list>self.runs[i]
                            chain = <Chain>run[j]
                            chain.advance1()
            else:
                r = self.mpiWorldRank
                i = (r % nchains) / self._ncoupled
                j = (r % nchains) % self._ncoupled
                run = <list>self.runs[i]
                chain = <Chain>run[j]
                chain.advance0()
                chain.advance1()

    cdef void advance(self) except *:
        IF @enable_mpi@:
            if self.mpiWorldSize == 1:
                self.advanceUni()
            else:
                self.advanceMpi()
        ELSE:
            self.advanceUni()

    # Compute Rcov convergence diagnostic.
    cdef double computeRcovUni(self, uint64_t step) except *:
        cdef ret
        cdef uint64_t first, past, lower, upper, i, j, n
        cdef unsigned runInd
        cdef double alphaEmp
        cdef double *rcovScratch = self.lnLs[self._nruns]
        cdef uint64_t last = step / self._stride

        # Utilize (at most) the last half of each chain.
        past = last + 1
        first = (last/2)+1
        if first == past:
            return 0.0

        # Compute the lower and upper bounds for the credibility intervals.  In
        # cases of rounding, expand rather than contract the credibility
        # intervals, in order to be conservative.  Expansion is conservative
        # because it is (for example) like computing the 97% credibility
        # interval instead of the 95% credibility interval, which is a more
        # stringent measure of agreement between sample sets.
        lower = <uint64_t>floor((self._cvgAlpha/2.0) * <double>(past-first))
        upper = (past-first-1) - lower

        ret = 0.0
        for 0 <= runInd < self._nruns:
            # Copy one run's samples into rcovScratch, and sort to make the
            # bounds of the credibility interval available.
            memcpy(rcovScratch, &self.lnLs[runInd][first], \
              (past-first)*sizeof(double))
            qsort(rcovScratch, past - first, sizeof(double), _lnLCmp)

            # Compute the proportion of the relevant lnLs that fall within the
            # credibility interval.
            n = 0
            for 0 <= i < self._nruns:
                for first <= j < past:
                    if self.lnLs[i][j] >= rcovScratch[lower] and \
                      self.lnLs[i][j] <= rcovScratch[upper]:
                        n += 1
            ret += <double>n / <double>(self._nruns * (past - first))
        ret /= <double>self._nruns
        # Compute the empirical alpha, based on the credibility interval
        # proportion, and use this to rescale Rcov, so that results can be
        # interpreted in the context of the intended alpha.
        alphaEmp = <double>(lower*2) / <double>(past-first)
        ret *= (1.0 - self._cvgAlpha) / (1.0 - alphaEmp)
        return ret

    IF @enable_mpi@:
        cdef double computeRcovMpi(self, uint64_t step) except *:
            cdef double rcov

            if self.mpiLeaderRank == 0:
                rcov = self.computeRcovUni(step)
            mpi.MPI_Bcast(&rcov, 1, mpi.MPI_DOUBLE, 0, mpi.MPI_COMM_WORLD)

            return rcov

    cdef double computeRcov(self, uint64_t step) except *:
        IF @enable_mpi@:
            if self.mpiWorldSize == 1:
                return self.computeRcovUni(step)
            else:
                return self.computeRcovMpi(step)
        ELSE:
            return self.computeRcovUni(step)

    cdef bint writeGraph(self, uint64_t step) except *:
        cdef file gfile
        cdef list cols0, cols1, lnLs0, lnLs1, colors
        cdef uint64_t samp, past, xsp, i
        cdef unsigned runInd
        cdef double lnL, rcov

        samp = step / self._stride
        past = samp+1
        if past < 4:
            # The code below does not generate an informative graph until there
            # are at least four samples to work with.
            return True

        # Compute which sample will be at the edge of the right graph.
        xsp = (samp/2)+1
        assert xsp > 0
        assert xsp < self.lnLsMax

        # Write to a temporary file, in order to be able to make the complete
        # file appear atomically in its final location.
        gfile = open("%s.lnL.R.tmp" % self.outPrefix, "w")

        cols0 = []
        cols1 = []
        for 0 <= runInd < self._nruns:
            lnLs0 = []
            lnLs1 = []
            for 0 <= i < xsp:
                lnLs0.append(self.lnLs[runInd][i])
            for xsp <= i < past:
                lnLs1.append(self.lnLs[runInd][i])
            cols0.append(", ".join(["%.6e" % lnL for lnL in lnLs0]))
            cols1.append(", ".join(["%.6e" % lnL for lnL in lnLs1]))
        gfile.write("lnLs0 = matrix(c(%s), ncol=%d)\n" % \
          (",\n".join(cols0), self._nruns))
        gfile.write("lnLs1 = matrix(c(%s), ncol=%d)\n" % \
          (",\n".join(cols1), self._nruns))

        gfile.write("x0 = seq(0, %d, by=%d)\n" % \
          (self._stride * (xsp-1), self._stride))
        gfile.write("x1 = seq(%d, %d, by=%d)\n" % \
          (self._stride * xsp, self._stride * samp, self._stride))

        gfile.write("par(mfrow=c(1, 2))\n")

        colors = ["black", "red", "green", "blue", "orange", "brown", "cyan", \
          "purple"]

        gfile.write('plot(c(min(x0), max(x0)), c(min(lnLs0), max(lnLs0)), pch="", xlab="step", ylab="lnL")\n')
        for 0 <= i < self._nruns:
            gfile.write('lines(x0, lnLs0[,%d], col="%s")\n' % \
              (i+1, colors[i % len(colors)]))

        gfile.write('plot(c(min(x1), max(x1)), c(min(lnLs1), max(lnLs1)), pch="", xlab="step", ylab="")\n')
        for 0 <= i < self._nruns:
            gfile.write('lines(x1, lnLs1[,%d], col="%s")\n' % \
              (i+1, colors[i % len(colors)]))

        gfile.close()
        os.rename("%s.lnL.R.tmp" % self.outPrefix, "%s.lnL.R" % self.outPrefix)

        return False

    cdef str formatRclass(self, Lik lik, unsigned model):
        cdef list strs, rclass
        cdef unsigned rlen, i

        strs = []
        strs.append("( ")
        rclass = lik.getRclass(model)
        rlen = len(rclass)
        for 0 <= i < rlen:
            if i > 0 and rlen > 10:
                strs.append(",")
            strs.append("%d" % rclass[i])
        strs.append(" )")
        return "".join(strs)

    cdef str formatRates(self, Lik lik, unsigned model, str fmt):
        cdef list strs, rclass
        cdef unsigned nrates, rlen, i
        cdef double r

        strs = []
        strs.append("[ ")
        nrates = lik.getNrates(model)
        rclass = lik.getRclass(model)
        rlen = len(rclass)
        for 0 <= i < rlen:
            if i > 0:
                strs.append(" ")
            r = lik.getRate(model, <unsigned>rclass[i])
            strs.append(fmt % r)
        strs.append(" ]")
        return "".join(strs)

    cdef str formatFreqs(self, Lik lik, unsigned model, str fmt):
        cdef list strs
        cdef unsigned nstates, i
        cdef double fSum, f

        strs = []
        strs.append("[ ")

        nstates = lik.char_.nstates
        fSum = 0.0
        for 0 <= i < nstates:
            f = lik.getFreq(model, i)
            fSum += f
        for 0 <= i < nstates:
            if i > 0:
                strs.append(" ")
            f = lik.getFreq(model, i) / fSum
            strs.append(fmt % f)

        strs.append(" ]")
        return "".join(strs)

    cdef str formatLnLs(self, uint64_t step, str fmt):
        cdef uint64_t samp
        cdef list strs
        cdef unsigned i

        samp = step / self._stride

        strs = []
        strs.append("[ ")
        for 0 <= i < self._nruns:
            if i > 0:
                strs.append(" ")
            strs.append(fmt % self.lnLs[i][samp])
        strs.append(" ]")
        return "".join(strs)

    cdef str formatRateStats(self, Mc3RateStats *rateStats):
        cdef list strs
        cdef unsigned i

        strs = []
        strs.append("[ ")
        for 0 <= i < self._nruns:
            if i > 0:
                strs.append(" ")
            strs.append("%.6f" % rateStats[i].ema)
        strs.append(" ]")
        return "".join(strs)

    cdef str formatPropStats(self):
        cdef list strs
        cdef unsigned i, j

        strs = []
        strs.append("{ ")
        for 0 <= i < PropCnt:
            if i > 0:
                strs.append(" ")
            if self.props[i] != 0.0:
                strs.append(self.formatRateStats(self.propStats[i]))
            else:
                strs.append("[ ")
                for 0 <= j < self._nruns:
                    if j > 0:
                        strs.append(" ")
                    strs.append("---")
                strs.append(" ]")
        strs.append(" }")

        return "".join(strs)

    # Write to .l log file.
    cdef void lWrite(self, str s) except *:
        IF @enable_mpi@:
            if self.mpiLeaderRank != 0:
                return

        self.lFile.write(s)
        self.lFile.flush()
        if self.verbose:
            sys.stdout.write(s)

    # Write to .t log file.
    cdef void tWrite(self, uint64_t step) except *:
        cdef unsigned i
        cdef Lik lik

        IF @enable_mpi@:
            if self.mpiLeaderRank != 0:
                return

        for 0 <= i < self._nruns:
            lik = <Lik>self.liks[i]
            self.tFile.write("[%d %d] " % (i, step))
            lik.tree.render(True, "%.12e", None, self.tFile)
        self.tFile.flush()

    # Write to .p log file.
    cdef void pWrite(self, uint64_t step) except *:
        cdef unsigned i, m
        cdef double wsum, w, wVar, wInvar, pinvar
        cdef Lik lik

        IF @enable_mpi@:
            if self.mpiLeaderRank != 0:
                return

        for 0 <= i < self._nruns:
            lik = <Lik>self.liks[i]
            wsum = 0.0
            nmodels = lik.nmodels()
            for 0 <= m < nmodels:
                wsum += lik.getWeight(m)
            assert wsum > 0.0
            for 0 <= m < nmodels:
                w = lik.getWeight(m)/wsum
                wVar = lik.getWVar(m)
                wInvar = lik.getWInvar(m)
                pinvar = wInvar / (wVar+wInvar)
                self.pFile.write( \
                  "%d\t%d\t%d\t%.5f\t%.5f\t%s%s %.5e %.5e %.5f %s\n" \
                  % (i, step, m, w, lik.getRmult(m), \
                  self.formatRclass(lik, m), self.formatRates(lik, m, "%.5e"), \
                  lik.getWNorm(), lik.getAlpha(m), pinvar, \
                  self.formatFreqs(lik, m, "%.5e")))
                if self.verbose:
                    sys.stdout.write( \
                      "p\t%d\t%d\t%d\t%.5f\t%.5f\t%s%s %.5f %.5e %.5f %s\n" % \
                      (i, step, m, w, lik.getRmult(m), \
                      self.formatRclass(lik, m), \
                      self.formatRates(lik, m, "%.4e"), lik.getWNorm(), \
                      lik.getAlpha(m), pinvar, \
                      self.formatFreqs(lik, m, "%.5f")))
        self.pFile.flush()

    # Write to .s log file.
    cdef void sWrite(self, uint64_t step, double rcov) \
      except *:
        cdef str swapStats, propStats, rcovStr

        IF @enable_mpi@:
            if self.mpiLeaderRank != 0:
                return

        rcovStr = ("%.6f" % rcov if rcov != -1.0 else "--------")
        swapStats = self.formatRateStats(self.swapStats)
        propStats = self.formatPropStats()
        self.sFile.write("%d\t%s %s %s %s\n" % (step, \
          self.formatLnLs(step, "%.11e"), rcovStr, swapStats, \
          propStats))
        self.sFile.flush()
        if self.verbose:
            sys.stdout.write("s\t%d\t%s %s\n" % (step, \
              self.formatLnLs(step, "%.6f"), rcovStr))

    cdef bint sample(self, uint64_t step) except *:
        cdef bint converged
        cdef double rcov
        cdef uint64_t samp = step / self._stride

        IF @enable_mpi@:
            self.storeSwapStats()
            self.storePropStats()
        self.updateDiags(step)
        self.storeLiksLnLs(step)

        if step < self._minStep or self._nruns == 1 or \
          step % (self._stride * self._cvgSampStride) != 0:
            rcov = -1.0
            converged = False
        else:
            rcov = self.computeRcov(step)
            converged = (rcov + self._cvgAlpha + self._cvgEpsilon >= 1.0)

        self.tWrite(step)
        self.pWrite(step)
        self.sWrite(step, rcov)

        self.resetLiks()

        return converged

    cdef void randomDnaQ(self, Lik lik, unsigned model, sfmt_t *prng) except *:
        """
            Initialize the specified model within lik's mixture vector with
            parameters drawn from their respective prior distributions.
        """
        cdef unsigned nstates, rlen, nrates, i
        cdef list rclasses, rclass
        cdef double u, norm, factor, cum
        global _dnaRclasses

        nstates = lik.char_.nstates
        rlen = nstates * (nstates-1) / 2

        if lik.nmodels() > 1:
            # Relative weight.
            if self.props[PropWeight] > 0.0:
                lik.setWeight(model, -log(1.0 - genrand_res53(prng)))

            # Rate multiplier.
            if self.props[PropRmult] > 0.0:
                lik.setRmult(model, -log(1.0 - genrand_res53(prng)))

        # State frequencies.
        if self.props[PropFreq] > 0.0:
            if self.props[PropFreqJump] > 0.0:
                u = genrand_res53(prng)
                # Equal/estimated frequencies, chosen according to the
                # frequency resolution prior (freqJumpPrior).
                norm = 0.0
                factor = 1.0
                for 0 <= i < 2:
                    norm += 1.0 / factor
                    factor *= self._freqJumpPrior
                cum = 0.0
                factor = 1.0
                for 1 <= i < 2:
                    cum += (1.0 / factor) / norm
                    if cum >= u:
                        break
                    factor *= self._freqJumpPrior
                if i == 2:
                    # Estimate frequencies.
                    for 0 <= i < nstates:
                        lik.setFreq(model, i, -log(1.0 - genrand_res53(prng)))
            else:
                for 0 <= i < nstates:
                    lik.setFreq(model, i, -log(1.0 - genrand_res53(prng)))

        # Relative mutation rates.
        if self.props[PropRate] > 0.0:
            if self.props[PropRateJump] > 0.0:
                u = genrand_res53(prng)
                # Number of rates, chosen according to the rclass resolution
                # prior (rateJumpPrior).
                norm = 0.0
                factor = 1.0
                for 0 <= i < rlen:
                    norm += 1.0 / factor
                    factor *= self._rateJumpPrior
                cum = 0.0
                factor = 1.0
                for 1 <= nrates < rlen+1:
                    cum += (1.0 / factor) / norm
                    if cum >= u:
                        break
                    factor *= self._rateJumpPrior
                # Rate class, randomly chosen from the appropriate resolution
                # class.
                if _dnaRclasses is None:
                    _initRclasses()
                rclasses = _dnaRclasses[nrates-1]
                rclass = rclasses[gen_rand64_range(prng, len(rclasses))]
            else:
                nrates = rlen
                rclass = range(rlen)
            lik.setRclass(model, rclass)
            for 0 <= i < nrates:
                lik.setRate(model, i, -log(1.0 - genrand_res53(prng)))

        if self.props[PropRateShapeInv] > 0.0:
            # Gamma-distributed rates shape parameter.
            if self._ncat > 1:
                if self.props[PropRateShapeInvJump] > 0.0:
                    u = genrand_res53(prng)
                    # Flat/+G, chosen according to the +G resolution prior
                    # (rateShapeInvJumpPrior).
                    norm = 0.0
                    factor = 1.0
                    for 0 <= i < 2:
                        norm += 1.0 / factor
                        factor *= self._rateShapeInvJumpPrior
                    cum = 0.0
                    factor = 1.0
                    for 1 <= i < 2:
                        cum += (1.0 / factor) / norm
                        if cum >= u:
                            break
                        factor *= self._rateShapeInvJumpPrior
                    if i == 2:
                        # +G.
                        lik.setAlpha(model, -log(1.0 - genrand_res53(prng)) \
                          * self._rateShapeInvPrior)
                else:
                    lik.setAlpha(model, -log(1.0 - genrand_res53(prng)) \
                      * self._rateShapeInvPrior)

        if self.props[PropInvar] > 0.0:
            # +I variable/invariable relative weight parameters.
            if self._invar:
                if self.props[PropInvarJump] > 0.0:
                    u = genrand_res53(prng)
                    # [not +I]/+I, chosen according to the +I resolution prior
                    # (invarJumpPrior).
                    norm = 0.0
                    factor = 1.0
                    for 0 <= i < 2:
                        norm += 1.0 / factor
                        factor *= self._invarJumpPrior
                    cum = 0.0
                    factor = 1.0
                    for 1 <= i < 2:
                        cum += (1.0 / factor) / norm
                        if cum >= u:
                            break
                        factor *= self._invarJumpPrior
                    if i == 2:
                        # +I.
                        lik.setWVar(model, -log(1.0 - genrand_res53(prng)) \
                          * (1.0 - self._invarPrior))
                        lik.setWInvar(model, -log(1.0 - genrand_res53(prng)) \
                          * self._invarPrior)
                else:
                    lik.setWVar(model, -log(1.0 - genrand_res53(prng)) \
                      * (1.0 - self._invarPrior))
                    lik.setWInvar(model, -log(1.0 - genrand_res53(prng)) \
                      * self._invarPrior)

    cpdef Lik randomLik(self, Tree tree=None):
        """
            Generate a Lik instance with all parameters drawn from their
            respective prior distributions.  By default, this method is used
            internally to generate independent starting points for Mc3 chains,
            but it can be used directly as the basis for manually creating
            starting points with some fixed parameters.  In order to fix
            parameters:

              1) Disable the appropriate proposal(s).
              2) Create a list of Lik instances via iterative calls to this
                 method.
              3) For each Lik instance, override the to-be-fixed parameter(s).
              4) Pass the list of Lik instances to the run() method.
        """
        cdef Lik lik
        cdef Edge edge
        cdef unsigned m, i
        cdef sfmt_t *prng

        prng = init_gen_rand(random.randint(0, 0xffffffff))
        if prng == NULL:
            raise MemoryError("Error allocating prng")

        try:
            # Generate a random fully resolved tree.  Polytomous starting trees
            # are never generated, but that is okay since there is no
            # requirement to draw the starting tree from the prior.  The main
            # goal here is to avoid systematic starting point dependence.
            if tree is None:
                tree = Tree(self.alignment.taxaMap.ntaxa, \
                  self.alignment.taxaMap)
                tree.deroot()

            # Branch lengths.
            for edge in tree.getEdges():
                edge.length = -log(1.0 - genrand_res53(prng)) / self._brlenPrior

            if self.props[PropMixtureJump] > 0.0:
                # Randomly draw _nmodels from the geometrically distributed
                # prior.
                self._nmodels = 1
                while genrand_res53(prng) > self._mixtureJumpPrior:
                    self._nmodels += 1

            lik = Lik(tree, self.alignment, self._nmodels, self._ncat, \
              self._catMedian, self._invar)

            # Randomly draw model parameters from their prior distributions.
            for 0 <= m < self._nmodels:
                self.randomDnaQ(lik, m, prng)
        finally:
            fini_gen_rand(prng)

        return lik

    cpdef run(self, bint verbose=False, list liks=None):
        """
            Run until convergence is reached, or the maximum number of steps
            is reached, whichever comes first.  Collate the results from the
            unheated chain(s) and write the raw results to disk.

            In order to specify initial parameter values, specify a list of Lik
            instances via 'liks'.
        """
        cdef double graphT0, graphT1
        cdef uint64_t step
        cdef bint graph

        self.verbose = verbose

        IF @enable_mpi@:
            self.initComms()
        try:
            self.initLogs()
            self.initSwapInfo()
            self.initSwapStats()
            self.initPropStats()
            self.initLiks()
            self.initLnLs()
            self.initProps()
            self.initRuns(liks)

            IF @enable_mpi@:
                if self.mpiLeaderRank == 0:
                    graph = (self._graphDelay >= 0.0)
                else:
                    graph = False
            ELSE:
                graph = (self._graphDelay >= 0.0)
            if graph:
                graphT0 = 0.0

            self.sample(0)
            # Run the chains no further than than _maxStep.
            for 1 <= step <= self._maxStep:
                self.advance()
                if step % self._stride == 0:
                    # sample() will not claim convergence until step is at
                    # least _minStep.
                    if self.sample(step):
                        self.lWrite("Runs converged\n")
                        break

                    if graph:
                        graphT1 = time.time()
                        if graphT1 - graphT0 >= self._graphDelay:
                            if not self.writeGraph(step):
                                graphT0 = time.time()

            # Write a graph one last time, regardless of whether graphs were
            # written previously.
            IF @enable_mpi@:
                if self.mpiLeaderRank == 0:
                    self.writeGraph(step)
            ELSE:
                self.writeGraph(step)
        except:
            error = sys.exc_info()
            IF @enable_mpi@:
                if self.mpiWorldSize > 1:
                    import traceback
                    sys.stderr.write("MPI node %s, rank %d of %d:" \
                      " Premature termination at %s: Exception %r\n" % \
                      (MPI.Get_processor_name(), MPI.COMM_WORLD.Get_rank(), \
                      MPI.COMM_WORLD.Get_size(), \
                      time.strftime("%Y/%m/%d %H:%M:%S (%Z)", \
                      time.localtime(time.time())), sys.exc_info()[1]))
                    for l in traceback.format_exception(*error):
                        sys.stderr.write(l)
                    mpi.MPI_Abort(mpi.MPI_COMM_WORLD, 1)
            ELSE:
                self.lWrite("Premature termination: Exception %r\n" % \
                  sys.exc_info()[1])
            raise
        finally:
            self.lWrite("Finish run: %s\n" % \
              time.strftime("%Y/%m/%d %H:%M:%S (%Z)", \
              time.localtime(time.time())))

    cdef double getGraphDelay(self):
        return self._graphDelay
    cdef void setGraphDelay(self, double graphDelay):
        self._graphDelay = graphDelay
    property graphDelay:
        """
            If non-negative, periodically (no more often than every graphDelay
            seconds) output an R program to <outPrefix>.lnL.R that generates a
            graphical representation of the log-likelihoods.  To monitor
            progress, use something like the following at an interactive R
            prompt:

              while (1) {source("outPrefix.lnL.R"); Sys.sleep(10)}

            A final graph will always be written, regardless of this setting.
        """
        def __get__(self):
            return self.getGraphDelay()
        def __set__(self, double graphDelay):
            self.setGraphDelay(graphDelay)

    cdef double getEmaAlpha(self):
        return self._emaAlpha
    cdef void setEmaAlpha(self, double emaAlpha) except *:
        if not (0.0 < emaAlpha and emaAlpha <= 1.0):
            raise ValueError("Validation failure: 0.0 < emaAlpha <= 1.0")
        self._emaAlpha = emaAlpha
    property emaAlpha:
        """
            Various diagnostic statistics (namely heat swap rates and proposal
            acceptance rates) are reported as exponential moving averages,
            using the following formulas:

              s(0) = x(0)
              s(t) = [emaAlpha * x(t)] + [s(t-1) * (1-emaAlpha)]
                   = s(t-1) + emaAlpha * [x(t) - s(t-1)]

            Discrete time (t) is measured in samples (one per stride).
            (0 < emaAlpha <= 1).  Larger values cause faster decay.
        """
        def __get__(self):
            return self.getEmaAlpha()
        def __set__(self, double emaAlpha):
            self.setEmaAlpha(emaAlpha)

    cdef uint64_t getCvgSampStride(self):
        return self._cvgSampStride
    cdef void setCvgSampStride(self, uint64_t cvgSampStride):
        if not (cvgSampStride >= 1):
            raise ValueError("Validation failure: cvgSampStride >= 1")
        self._cvgSampStride = cvgSampStride
    property cvgSampStride:
        """
            cvgSampStride limits convergence diagnostics computation to once
            every stride*cvgSampStride steps.  This is primarily useful for
            reducing convergence diagnostic overhead.
        """
        def __get__(self):
            return self.getCvgSampStride()
        def __set__(self, uint64_t cvgSampStride):
            self.setCvgSampStride(cvgSampStride)

    cdef double getCvgAlpha(self):
        return self._cvgAlpha
    cdef void setCvgAlpha(self, double cvgAlpha) except *:
        if not (0.0 <= cvgAlpha and cvgAlpha <= 1.0):
            raise ValueError("Validation failure: 0.0 <= cvgAlpha <= 1.0")
        self._cvgAlpha = cvgAlpha
    property cvgAlpha:
        """
            (1.0 - cvgAlpha) is the confidence interval used for diagnosing
            convergence.
        """
        def __get__(self):
            return self.getCvgAlpha()
        def __set__(self, double cvgAlpha):
            self.setCvgAlpha(cvgAlpha)

    cdef double getCvgEpsilon(self):
        return self._cvgEpsilon
    cdef void setCvgEpsilon(self, double cvgEpsilon) except *:
        if not (0.0 <= self._cvgEpsilon and self._cvgEpsilon <= 1.0):
            raise ValueError("Validation failure: 0.0 <= cvgEpsilon <= 1.0")
        self._cvgEpsilon = cvgEpsilon
    property cvgEpsilon:
        """
            Convergence is achieved once the mean coverage is within cvgEpsilon
            of the confidence interval.
        """
        def __get__(self):
            return self.getCvgEpsilon()
        def __set__(self, double cvgEpsilon):
            self.setCvgEpsilon(cvgEpsilon)

    cdef uint64_t getMinStep(self):
        return self._minStep
    cdef void setMinStep(self, uint64_t minStep) except *:
        if not (minStep <= self._maxStep):
            raise ValueError("Validation failure: minStep <= maxStep")
        self._minStep = minStep
    property minStep:
        """
            Minimum final step index for each chain.
        """
        def __get__(self):
            return self.getMinStep()
        def __set__(self, uint64_t minStep):
            self.setMinStep(minStep)

    cdef uint64_t getMaxStep(self):
        return self._maxStep
    cdef void setMaxStep(self, uint64_t maxStep) except *:
        if not (self._minStep <= maxStep):
            raise ValueError("Validation failure: minStep <= maxStep")
        self._maxStep = maxStep
    property maxStep:
        """
            Maximum final step index for each chain.
        """
        def __get__(self):
            return self.getMaxStep()
        def __set__(self, uint64_t maxStep):
            self.setMaxStep(maxStep)

    cdef unsigned getStride(self):
        return self._stride
    cdef void setStride(self, unsigned stride) except *:
        if not stride >= 1:
            raise ValueError("Validation failure: stride >= 1")
        self._stride = stride
    property stride:
        """
            Chain sample interval.
        """
        def __get__(self):
            return self.getStride()
        def __set__(self, unsigned stride):
            self.setStride(stride)

    cdef unsigned getNruns(self):
        return self._nruns
    cdef void setNruns(self, unsigned nruns) except *:
        if not nruns >= 1:
            raise ValueError("Validation failure: nruns >= 1")
        self._nruns = nruns
    property nruns:
        """
            Number of independent runs.  If there is only one run, convergence
            diagnostics will not be computed.
        """
        def __get__(self):
            return self.getNruns()
        def __set__(self, unsigned nruns):
            self.setNruns(nruns)

    cdef unsigned getNcoupled(self):
        return self._ncoupled
    cdef void setNcoupled(self, unsigned ncoupled) except *:
        if not ncoupled >= 1:
            raise ValueError("Validation failure: ncoupled >= 1")
        self._ncoupled = ncoupled
    property ncoupled:
        """
            Number of Metropolis-coupled chains per run.
        """
        def __get__(self):
            return self.getNcoupled()
        def __set__(self, unsigned ncoupled):
            self.setNcoupled(ncoupled)

    cdef double getHeatDelta(self):
        return self._heatDelta
    cdef void setHeatDelta(self, double heatDelta) except *:
        if not heatDelta > 0.0:
            raise ValueError("Validation failure: heatDelta > 0.0")
        self._heatDelta = heatDelta
    property heatDelta:
        """
            Temperature interval between Metropolis-coupled chains.
        """
        def __get__(self):
            return self.getHeatDelta()
        def __set__(self, double heatDelta):
            self.setHeatDelta(heatDelta)

    cdef unsigned getSwapStride(self):
        return self._swapStride
    cdef void setSwapStride(self, unsigned swapStride) except *:
        if not swapStride >= 1:
            raise ValueError("Validation failure: swapStride >= 1")
        self._swapStride = swapStride
    property swapStride:
        """
            Chain heat swap interval.
        """
        def __get__(self):
            return self.getSwapStride()
        def __set__(self, unsigned swapStride):
            self.setSwapStride(swapStride)

    cdef unsigned getNmodels(self):
        return self._nmodels
    cdef void setNmodels(self, unsigned nmodels) except *:
        if not nmodels >= 1:
            raise ValueError("Validation failure: nmodels >= 1")
        self._nmodels = nmodels
    property nmodels:
        """
            Number of mixture models.  This setting is irrelevant if mixture
            jumps are enabled (default).
        """
        def __get__(self):
            return self.getNmodels()
        def __set__(self, unsigned nmodels):
            self.setNmodels(nmodels)

    cdef unsigned getNcat(self):
        return self._ncat
    cdef void setNcat(self, unsigned ncat) except *:
        if not ncat >= 1:
            raise ValueError("Validation failure: ncat >= 1")
        self._ncat = ncat
    property ncat:
        """
            Number of discrete Gamma-distributed relative mutation rate
            categories.
        """
        def __get__(self):
            return self.getNcat()
        def __set__(self, unsigned ncat):
            self.setNcat(ncat)

    cdef bint getCatMedian(self):
        return self._catMedian
    cdef void setCatMedian(self, bint catMedian):
        self._catMedian = catMedian
    property catMedian:
        """
            Use category means for Gamma-distributed relative mutation rate
            discretization if false; use category medians otherwise.
        """
        def __get__(self):
            return self.getCatMedian()
        def __set__(self, bint catMedian):
            self.setCatMedian(catMedian)

    cdef bint getInvar(self):
        return self._invar
    cdef void setInvar(self, bint invar):
        self._invar = invar
    property invar:
        """
            Support a proportion of invariable sites if true.
        """
        def __get__(self):
            return self.getInvar()
        def __set__(self, bint invar):
            self.setInvar(invar)

    cdef double getWeightLambda(self):
        return self._weightLambda
    cdef void setWeightLambda(self, double weightLambda) except *:
        if not weightLambda >= 0.0:
            raise ValueError("Validation failure: weightLambda >= 0.0")
        self._weightLambda = weightLambda
    property weightLambda:
        """
            Relative model weight multiplier.  This controls the range of
            weight change for weight change proposals.

            Proposed weights are drawn from a log-transformed sliding
            window.  A multiplier, m, is drawn from

                         1
            g(m) = ----------------
                   weightLambda * m

                            /        1             weightLambda/2 \ 
            in the interval | ----------------- , e               | .
                            |   weightLambda/2                    |
                            \  e                                  /

            Where weightLambda = 2 * log(a), proposed weights w* are in
            the interval (w/a, aw).
        """
        def __get__(self):
            return self.getWeightLambda()
        def __set__(self, double weightLambda):
            self.setWeightLambda(weightLambda)

    cdef double getFreqLambda(self):
        return self._freqLambda
    cdef void setFreqLambda(self, double freqLambda) except *:
        if not freqLambda >= 0.0:
            raise ValueError("Validation failure: freqLambda >= 0.0")
        self._freqLambda = freqLambda
    property freqLambda:
        """
            State frequency multiplier.  This controls the range of state
            frequency change for state frequency change proposals.

            Proposed state frequencies are drawn from a log-transformed sliding
            window.  A multiplier, m, is drawn from

                         1
            g(m) = --------------
                   freqLambda * m

                            /        1           freqLambda/2 \ 
            in the interval | --------------- , e             | .
                            |   freqLambda/2                  |
                            \  e                              /

            Where freqLambda = 2 * log(a), proposed frequencies f* are in
            the interval (f/a, af).
        """
        def __get__(self):
            return self.getFreqLambda()
        def __set__(self, double freqLambda):
            self.setFreqLambda(freqLambda)

    cdef double getRmultLambda(self):
        return self._rmultLambda
    cdef void setRmultLambda(self, double rmultLambda) except *:
        if not rmultLambda >= 0.0:
            raise ValueError("Validation failure: rmultLambda >= 0.0")
        self._rmultLambda = rmultLambda
    property rmultLambda:
        """
            Rate multiplier multiplier (yes, really).  This controls the range
            of rate multiplier change for rate multiplier change proposals.

            Proposed rate multipliers are drawn from a log-transformed sliding
            window.  A multiplier, m, is drawn from

                         1
            g(m) = ---------------
                   rmultLambda * m

                            /        1           rmultLambda/2 \ 
            in the interval | --------------- , e              | .
                            |   rmultLambda/2                  |
                            \  e                               /

            Where rmultLambda = 2 * log(a), proposed rate multipliers r* are in
            the interval (r/a, ar).
        """
        def __get__(self):
            return self.getRmultLambda()
        def __set__(self, double rmultLambda):
            self.setRmultLambda(rmultLambda)

    cdef double getRateLambda(self):
        return self._rateLambda
    cdef void setRateLambda(self, double rateLambda) except *:
        if not rateLambda >= 0.0:
            raise ValueError("Validation failure: rateLambda >= 0.0")
        self._rateLambda = rateLambda
    property rateLambda:
        """
            Relative mutation rate multiplier.  This controls the range of
            rate change for rate change proposals.

            Proposed rates are drawn from a log-transformed sliding window.  A
            multiplier, m, is drawn from

                         1
            g(m) = ---------------
                   rateLambda * m

                            /        1           rateLambda/2 \ 
            in the interval | --------------- , e             | .
                            |   rateLambda/2                  |
                            \  e                              /

            Where rateLambda = 2 * log(a), proposed rates r* are in the
            interval (r/a, ar).
        """
        def __get__(self):
            return self.getRateLambda()
        def __set__(self, double rateLambda):
            self.setRateLambda(rateLambda)

    cdef double getRateShapeInvLambda(self):
        return self._rateShapeInvLambda
    cdef void setRateShapeInvLambda(self, double rateShapeInvLambda) except *:
        if not rateShapeInvLambda >= 0.0:
            raise ValueError("Validation failure: rateShapeInvLambda >= 0.0")
        self._rateShapeInvLambda = rateShapeInvLambda
    property rateShapeInvLambda:
        """
            Gamma-distributed inverse shape change multiplier.  This controls
            the range of inverse shape change for inverse shape change
            proposals.

            Proposed inverse shapes are drawn from a log-transformed sliding
            window.  A multiplier, m, is drawn from

                             1
            g(m) = ----------------------
                   rateShapeInvLambda * m

                            /        1                  rateShapeInvLambda/2 \ 
            in the interval | ---------------------- , e                     | .
                            |   rateShapeInvLambda/2                         |
                            \  e                                             /

            Where rateShapeInvLambda = 2 * log(a), proposed inverse shapes s*
            are in the interval (s/a, as).
        """
        def __get__(self):
            return self.getRateShapeInvLambda()
        def __set__(self, double rateShapeInvLambda):
            self.setRateShapeInvLambda(rateShapeInvLambda)

    cdef double getInvarLambda(self):
        return self._invarLambda
    cdef void setInvarLambda(self, double invarLambda) except *:
        if not invarLambda >= 0.0:
            raise ValueError("Validation failure: invarLambda >= 0.0")
        self._invarLambda = invarLambda
    property invarLambda:
        """
            Variable/invariable site relative weight multiplier.  This controls
            the range of weight change for variable/invariable weight change
            proposals.  Variable and invariable site weights are exponentially
            distributed, which is equivalent in effect to using a flat
            Dirichlet prior for the proportion of invariant sites.

            Proposed weights are drawn from a log-transformed sliding
            window.  A multiplier, m, is drawn from

                         1
            g(m) = ---------------
                   invarLambda * m

                            /        1             invarLambda/2 \ 
            in the interval | ----------------- , e              | .
                            |   invarLambda/2                    |
                            \  e                                 /

            Where invarLambda = 2 * log(a), proposed weights w* are in the
            interval (w/a, aw).
        """
        def __get__(self):
            return self.getInvarLambda()
        def __set__(self, double invarLambda):
            self.setInvarLambda(invarLambda)

    cdef double getBrlenLambda(self):
        return self._brlenLambda
    cdef void setBrlenLambda(self, double brlenLambda) except *:
        if not brlenLambda >= 0.0:
            raise ValueError("Validation failure: brlenLambda >= 0.0")
        self._brlenLambda = brlenLambda
    property brlenLambda:
        """
            Branch change length multiplier.  This controls the range of branch
            length change for basic branch length change proposals.

            Proposed branch lengths are drawn from a log-transformed sliding
            window.  A multiplier, m, is drawn from

                         1
            g(m) = ---------------
                   brlenLambda * m

                            /        1            brlenLambda/2 \ 
            in the interval | ---------------- , e              | .
                            |   brlenLambda/2                   |
                            \  e                                /

            Where brlenLambda = 2 * log(a), proposed branch lengths v* are in
            the interval (v/a, av).
        """
        def __get__(self):
            return self.getBrlenLambda()
        def __set__(self, double brlenLambda):
            self.setBrlenLambda(brlenLambda)

    cdef double getEtbrPExt(self):
        return self._etbrPExt
    cdef void setEtbrPExt(self, double etbrPExt) except *:
        if not (0.0 <= etbrPExt and etbrPExt <= 1.0):
            raise ValueError("Validation failure: 0.0 <= etbrPExt <= 1.0")
        self._etbrPExt = etbrPExt
    property etbrPExt:
        """
            eTBR extension probability.  This controls how far away from the
            bisection edge reattachment occurs on average.
        """
        def __get__(self):
            return self.getEtbrPExt()
        def __set__(self, double etbrPExt):
            self.setEtbrPExt(etbrPExt)

    cdef double getEtbrLambda(self):
        return self._etbrLambda
    cdef void setEtbrLambda(self, double etbrLambda) except *:
        if not etbrLambda >= 0.0:
            raise ValueError("Validation failure: etbrLambda >= 0.0")
        self._etbrLambda = etbrLambda
    property etbrLambda:
        """
            eTBR branch length multiplier.  This controls the range of branch
            length change.
        """
        def __get__(self):
            return self.getEtbrLambda()
        def __set__(self, double etbrLambda):
            self.setEtbrLambda(etbrLambda)

    cdef double getRateShapeInvPrior(self):
        return self._rateShapeInvPrior
    cdef void setRateShapeInvPrior(self, double rateShapeInvPrior) except *:
        if not rateShapeInvPrior > 0.0:
            raise ValueError("Validation failure: rateShapeInvPrior > 0.0")
        self._rateShapeInvPrior = rateShapeInvPrior
    property rateShapeInvPrior:
        """
            Prior for exponentially distributed inverse (variance) of shape
            parameter for Gamma-distributed relative mutation rates.
        """
        def __get__(self):
            return self.getRateShapeInvPrior()
        def __set__(self, double rateShapeInvPrior):
            self.setRateShapeInvPrior(rateShapeInvPrior)

    cdef double getInvarPrior(self):
        return self._invarPrior
    cdef void setInvarPrior(self, double invarPrior) except *:
        if not invarPrior > 0.0:
            raise ValueError("Validation failure: invarPrior > 0.0")
        self._invarPrior = invarPrior
    property invarPrior:
        """
            Prior mean for proportion of invariable sites.
        """
        def __get__(self):
            return self.getInvarPrior()
        def __set__(self, double invarPrior):
            self.setInvarPrior(invarPrior)

    cdef double getBrlenPrior(self):
        return self._brlenPrior
    cdef void setBrlenPrior(self, double brlenPrior) except *:
        if not brlenPrior > 0.0:
            raise ValueError("Validation failure: brlenPrior > 0.0")
        self._brlenPrior = brlenPrior
    property brlenPrior:
        """
            Prior for exponentially distributed branch lengths.
        """
        def __get__(self):
            return self.getBrlenPrior()
        def __set__(self, double brlenPrior):
            self.setBrlenPrior(brlenPrior)

    cdef double getRateJumpPrior(self):
        return self._rateJumpPrior
    cdef void setRateJumpPrior(self, double rateJumpPrior) except *:
        if not rateJumpPrior > 0.0:
            raise ValueError("Validation failure: rateJumpPrior > 0.0")
        self._rateJumpPrior = rateJumpPrior
    property rateJumpPrior:
        """
            Prior for number of relative mutation rate parameters.  1.0
            indicates a flat prior; <1.0 favors more complex models, and >1.0
            favors less complex models.
        """
        def __get__(self):
            return self.getRateJumpPrior()
        def __set__(self, double rateJumpPrior):
            self.setRateJumpPrior(rateJumpPrior)

    cdef double getPolytomyJumpPrior(self):
        return self._polytomyJumpPrior
    cdef void setPolytomyJumpPrior(self, double polytomyJumpPrior) except *:
        if not polytomyJumpPrior > 0.0:
            raise ValueError("Validation failure: polytomyJumpPrior > 0.0")
        self._polytomyJumpPrior = polytomyJumpPrior
    property polytomyJumpPrior:
        """
            Prior for number of internal branches.  1.0 indicates a flat prior
            (all tree topologies, polytomous and fully resolved are given equal
            weight), <1.0 favors more resolved trees, and >1.0 favors less
            resolved trees.
        """
        def __get__(self):
            return self.getPolytomyJumpPrior()
        def __set__(self, double polytomyJumpPrior):
            self.setPolytomyJumpPrior(polytomyJumpPrior)

    cdef double getRateShapeInvJumpPrior(self):
        return self._rateShapeInvJumpPrior
    cdef void setRateShapeInvJumpPrior(self, double rateShapeInvJumpPrior) \
      except *:
        if not rateShapeInvJumpPrior > 0.0:
            raise ValueError("Validation failure: rateShapeInvJumpPrior > 0.0")
        self._rateShapeInvJumpPrior = rateShapeInvJumpPrior
    property rateShapeInvJumpPrior:
        """
            Prior for presence/absence of the Gamma-distributed relative
            mutation rates shape parameter.  Such models are commonly referred
            to as +G (plus Gamma) models (e.g. GTR+G).  1.0 indicates a flat
            prior (no preference for/against +G models), <1.0 favors +G models,
            and >1.0 favors +G-less models.

            The rateShapeInvJumpPrior parameterization is chosen to be
            consistent with priors such as rateJumpPrior and polytomyJumpPrior.
            Conceptually, all three are structured the same way, but
            rateShapeInvJumpPrior only applies to two levels of nested models,
            whereas the other priors may apply to many levels.
        """
        def __get__(self):
            return self.getRateShapeInvJumpPrior()
        def __set__(self, double rateShapeInvJumpPrior):
            self.setRateShapeInvJumpPrior(rateShapeInvJumpPrior)

    cdef double getInvarJumpPrior(self):
        return self._invarJumpPrior
    cdef void setInvarJumpPrior(self, double invarJumpPrior) \
      except *:
        if not invarJumpPrior > 0.0:
            raise ValueError("Validation failure: invarJumpPrior > 0.0")
        self._invarJumpPrior = invarJumpPrior
    property invarJumpPrior:
        """
            Prior for presence/absence of the proportion of invariable sites
            parameter.  Such models are commonly referred to as +I models (e.g.
            GTR+I).  1.0 indicates a flat prior (no preference for/against +I
            models), <1.0 favors +I models, and >1.0 favors +I-less models.

            The invarJumpPrior parameterization is chosen to be consistent with
            priors such as rateJumpPrior and polytomyJumpPrior.  Conceptually,
            all three are structured the same way, but invarJumpPrior only
            applies to two levels of nested models, whereas the other priors
            may apply to many levels.
        """
        def __get__(self):
            return self.getInvarJumpPrior()
        def __set__(self, double invarJumpPrior):
            self.setInvarJumpPrior(invarJumpPrior)

    cdef double getFreqJumpPrior(self):
        return self._freqJumpPrior
    cdef void setFreqJumpPrior(self, double freqJumpPrior) except *:
        if not freqJumpPrior > 0.0:
            raise ValueError("Validation failure: freqJumpPrior > 0.0")
        self._freqJumpPrior = freqJumpPrior
    property freqJumpPrior:
        """
            Prior for estimated/equal state frequencies.  1.0 indicates a flat
            prior (no preference for/against estimated frequencies), <1.0
            favors estimated frequencies, and >1.0 favors equal frequencies.

            The freqJumpPrior parameterization is chosen to be consistent with
            priors such as rateJumpPrior and polytomyJumpPrior.  Conceptually,
            all three are structured the same way, but freqJumpPrior only
            applies to two levels of nested models, whereas the other priors
            may apply to many levels.
        """
        def __get__(self):
            return self.getFreqJumpPrior()
        def __set__(self, double freqJumpPrior):
            self.setFreqJumpPrior(freqJumpPrior)

    cdef double getMixtureJumpPrior(self):
        return self._mixtureJumpPrior
    cdef void setMixtureJumpPrior(self, double mixtureJumpPrior) except *:
        if not (0.0 <= mixtureJumpPrior and mixtureJumpPrior < 1.0):
            raise ValueError( \
              "Validation failure: 0.0 <= mixtureJumpPrior < 1.0")
        self._mixtureJumpPrior = mixtureJumpPrior
    property mixtureJumpPrior:
        """
            Prior inverse mean number of models in the mixture, assuming a
            geometric distribution.  The prior probability of mixture degree
            d(M) is defined as:
                                                k
              P(d(M)=k) = (1.0-mixtureJumpPrior) * mixtureJumpPrior
        """
        def __get__(self):
            return self.getMixtureJumpPrior()
        def __set__(self, double mixtureJumpPrior):
            self.setMixtureJumpPrior(mixtureJumpPrior)

    cdef double getWeightProp(self):
        return self.props[PropWeight]
    cdef void setWeightProp(self, double weightProp) except *:
        if not weightProp >= 0.0:
            raise ValueError("Validation failure: weightProp >= 0.0")
        self.props[PropWeight] = weightProp
    property weightProp:
        """
            Relative proportion of proposals that modify relative model weights.
        """
        def __get__(self):
            return self.getWeightProp()
        def __set__(self, double weightProp):
            self.setWeightProp(weightProp)

    cdef double getFreqProp(self):
        return self.props[PropFreq]
    cdef void setFreqProp(self, double freqProp) except *:
        if not freqProp >= 0.0:
            raise ValueError("Validation failure: freqProp >= 0.0")
        self.props[PropFreq] = freqProp
    property freqProp:
        """
            Relative proportion of proposals that modify state frequencies.
        """
        def __get__(self):
            return self.getFreqProp()
        def __set__(self, double freqProp):
            self.setFreqProp(freqProp)

    cdef double getRmultProp(self):
        return self.props[PropRmult]
    cdef void setRmultProp(self, double rmultProp) except *:
        if not rmultProp >= 0.0:
            raise ValueError("Validation failure: rmultProp >= 0.0")
        self.props[PropRmult] = rmultProp
    property rmultProp:
        """
            Relative proportion of proposals that modify rate multipliers.
        """
        def __get__(self):
            return self.getRmultProp()
        def __set__(self, double rmultProp):
            self.setRmultProp(rmultProp)

    cdef double getRateProp(self):
        return self.props[PropRate]
    cdef void setRateProp(self, double rateProp) except *:
        if not rateProp >= 0.0:
            raise ValueError("Validation failure: rateProp >= 0.0")
        self.props[PropRate] = rateProp
    property rateProp:
        """
            Relative proportion of proposals that modify relative mutation
            rates.
        """
        def __get__(self):
            return self.getRateProp()
        def __set__(self, double rateProp):
            self.setRateProp(rateProp)

    cdef double getRateShapeInvProp(self):
        return self.props[PropRateShapeInv]
    cdef void setRateShapeInvProp(self, double rateShapeInvProp) except *:
        if not rateShapeInvProp >= 0.0:
            raise ValueError("Validation failure: rateShapeInvProp >= 0.0")
        self.props[PropRateShapeInv] = rateShapeInvProp
    property rateShapeInvProp:
        """
            Relative proportion of proposals that modify inverse (variance) of
            shape parameters for Gamma-distributed relative mutation rates.
        """
        def __get__(self):
            return self.getRateShapeInvProp()
        def __set__(self, double rateShapeInvProp):
            self.setRateShapeInvProp(rateShapeInvProp)

    cdef double getInvarProp(self):
        return self.props[PropInvar]
    cdef void setInvarProp(self, double invarProp) except *:
        if not invarProp >= 0.0:
            raise ValueError("Validation failure: invarProp >= 0.0")
        self.props[PropInvar] = invarProp
    property invarProp:
        """
            Relative proportion of proposals that modify the proportion of
            invariable sites.
        """
        def __get__(self):
            return self.getInvarProp()
        def __set__(self, double invarProp):
            self.setInvarProp(invarProp)

    cdef double getBrlenProp(self):
        return self.props[PropBrlen]
    cdef void setBrlenProp(self, double brlenProp) except *:
        if not brlenProp >= 0.0:
            raise ValueError("Validation failure: brlenProp >= 0.0")
        self.props[PropBrlen] = brlenProp
    property brlenProp:
        """
            Relative proportion of proposals that modify branch lengths.
        """
        def __get__(self):
            return self.getBrlenProp()
        def __set__(self, double brlenProp):
            self.setBrlenProp(brlenProp)

    cdef double getEtbrProp(self):
        return self.props[PropEtbr]
    cdef void setEtbrProp(self, double etbrProp) except *:
        if not etbrProp >= 0.0:
            raise ValueError("Validation failure: etbrProp >= 0.0")
        self.props[PropEtbr] = etbrProp
    property etbrProp:
        """
            Relative proportion of proposals that modify tree topology via
            eTBR.
        """
        def __get__(self):
            return self.getEtbrProp()
        def __set__(self, double etbrProp):
            self.setEtbrProp(etbrProp)

    cdef double getRateJumpProp(self):
        return self.props[PropRateJump]
    cdef void setRateJumpProp(self, double rateJumpProp) except *:
        if not rateJumpProp >= 0.0:
            raise ValueError("Validation failure: rateJumpProp >= 0.0")
        self.props[PropRateJump] = rateJumpProp
    property rateJumpProp:
        """
            Relative proportion of proposals that split/merge relative mutation
            rate parameters.
        """
        def __get__(self):
            return self.getRateJumpProp()
        def __set__(self, double rateJumpProp):
            self.setRateJumpProp(rateJumpProp)

    cdef double getPolytomyJumpProp(self):
        return self.props[PropPolytomyJump]
    cdef void setPolytomyJumpProp(self, double polytomyJumpProp) except *:
        if not polytomyJumpProp >= 0.0:
            raise ValueError("Validation failure: polytomyJumpProp >= 0.0")
        self.props[PropPolytomyJump] = polytomyJumpProp
    property polytomyJumpProp:
        """
            Relative proportion of proposals that add/remove internal branches.
        """
        def __get__(self):
            return self.getPolytomyJumpProp()
        def __set__(self, double polytomyJumpProp):
            self.setPolytomyJumpProp(polytomyJumpProp)

    cdef double getRateShapeInvJumpProp(self):
        return self.props[PropRateShapeInvJump]
    cdef void setRateShapeInvJumpProp(self, double rateShapeInvJumpProp) \
      except *:
        if not rateShapeInvJumpProp >= 0.0:
            raise ValueError("Validation failure: rateShapeInvJumpProp >= 0.0")
        self.props[PropRateShapeInvJump] = rateShapeInvJumpProp
    property rateShapeInvJumpProp:
        """
            Relative proportion of proposals that add/remove Gamma-distributed
            ralative rates parameters (+G models).
        """
        def __get__(self):
            return self.getRateShapeInvJumpProp()
        def __set__(self, double rateShapeInvJumpProp):
            self.setRateShapeInvJumpProp(rateShapeInvJumpProp)

    cdef double getInvarJumpProp(self):
        return self.props[PropInvarJump]
    cdef void setInvarJumpProp(self, double invarJumpProp) \
      except *:
        if not invarJumpProp >= 0.0:
            raise ValueError("Validation failure: invarJumpProp >= 0.0")
        self.props[PropInvarJump] = invarJumpProp
    property invarJumpProp:
        """
            Relative proportion of proposals that add/remove invariable site
            proportion parameters (+I models).
        """
        def __get__(self):
            return self.getInvarJumpProp()
        def __set__(self, double invarJumpProp):
            self.setInvarJumpProp(invarJumpProp)

    cdef double getFreqJumpProp(self):
        return self.props[PropFreqJump]
    cdef void setFreqJumpProp(self, double freqJumpProp) except *:
        if not freqJumpProp >= 0.0:
            raise ValueError("Validation failure: freqJumpProp >= 0.0")
        self.props[PropFreqJump] = freqJumpProp
    property freqJumpProp:
        """
            Relative proportion of proposals that switch between
            estimated/equal state frequencies.
        """
        def __get__(self):
            return self.getFreqJumpProp()
        def __set__(self, double freqJumpProp):
            self.setFreqJumpProp(freqJumpProp)

    cdef double getMixtureJumpProp(self):
        return self.props[PropMixtureJump]
    cdef void setMixtureJumpProp(self, double mixtureJumpProp) except *:
        if not mixtureJumpProp >= 0.0:
            raise ValueError("Validation failure: mixtureJumpProp >= 0.0")
        self.props[PropMixtureJump] = mixtureJumpProp
    property mixtureJumpProp:
        """
            Relative proportion of proposals that add/remove Q matrices to/from
            the model mixture.
        """
        def __get__(self):
            return self.getMixtureJumpProp()
        def __set__(self, double mixtureJumpProp):
            self.setMixtureJumpProp(mixtureJumpProp)
