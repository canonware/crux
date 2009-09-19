# Forward declaration.
cdef class Chain

cdef enum:
    PropWeight           =  0
    PropFreq             =  1
    PropRmult            =  2
    PropRate             =  3
    PropRateShapeInv     =  4
    PropInvar            =  5
    PropBrlen            =  6
    PropEtbr             =  7
    PropRateJump         =  8
    PropPolytomyJump     =  9
    PropRateShapeInvJump = 10
    PropInvarJump        = 11
    PropFreqJump         = 12
    PropMixtureJump      = 13
    PropCnt              = 14 # Number of proposals.

from libc cimport uint64_t
from SFMT cimport sfmt_t
from Crux.Tree cimport Tree
from Crux.Tree.Bipart cimport Bipart
from Crux.Tree.Lik cimport Lik
from Crux.Mc3 cimport Mc3

cdef class Chain:
    cdef Mc3 master
    cdef unsigned run
    cdef unsigned ind
    cdef uint64_t nswap
    cdef uint64_t accepts[PropCnt]
    cdef uint64_t rejects[PropCnt]
    cdef double heat
    cdef unsigned swapInd
    cdef double swapProb
    cdef sfmt_t *swapPrng
    cdef sfmt_t *prng
    cdef Tree tree
    cdef Lik lik
    cdef double lnL
    cdef uint64_t step

    cdef bint weightPropose(self) except *
    cdef bint freqsEqual(self, Lik lik, unsigned mInd)
    cdef unsigned nModelsFreqsEstim(self, Lik lik)
    cdef bint freqPropose(self) except *
    cdef bint rmultPropose(self) except *
    cdef unsigned nModelsRatesEstim(self, Lik lik)
    cdef bint ratePropose(self) except *
    cdef unsigned nModelsRatesGamma(self, Lik lik)
    cdef bint rateShapeInvPropose(self) except *
    cdef unsigned nModelsInvar(self, Lik lik)
    cdef bint invarPropose(self) except *
    cdef bint brlenPropose(self) except *
    cdef bint etbrPropose(self) except *
    cdef void rateMergePropose(self, unsigned mInd, list rclass, \
      unsigned nrates) except *
    cdef void rateSplitPropose(self, unsigned mInd, list rclass, \
      unsigned nrates) except *
    cdef bint rateJumpPropose(self) except *
    cdef void polytomyMergePropose(self, Tree tree, unsigned nedges, \
      unsigned ntaxa) except *
    cdef void polytomySplitPropose(self, Tree tree, unsigned nedges, \
      unsigned ntaxa) except *
    cdef bint polytomyJumpPropose(self) except *
    cdef void rateShapeInvRemovePropose(self, unsigned mInd, double alpha0) \
      except *
    cdef void rateShapeInvAddPropose(self, unsigned mInd) except *
    cdef bint rateShapeInvJumpPropose(self) except *
    cdef void invarRemovePropose(self, unsigned mInd, double wVar, \
      double wInvar) except *
    cdef void invarAddPropose(self, unsigned mInd) except *
    cdef bint invarJumpPropose(self) except *
    cdef void freqEqualPropose(self, unsigned mInd) except *
    cdef void freqEstimPropose(self, unsigned mInd) except *
    cdef bint freqJumpPropose(self) except *
    cdef void mixtureRemovePropose(self, unsigned nmodels) except *
    cdef void mixtureAddPropose(self, unsigned nmodels) except *
    cdef bint mixtureJumpPropose(self) except *
    cdef void advance0(self) except *
    cdef void advance1(self) except *
