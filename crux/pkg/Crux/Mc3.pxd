from libc cimport uint64_t
from SFMT cimport sfmt_t
from Crux.CTMatrix cimport Alignment
from Crux.Tree cimport Tree
from Crux.Tree.Lik cimport Lik

# Forward declarations.
cdef class Mc3Chain
cdef class Mc3

cdef struct Mc3SwapInfo:
    uint64_t step
    double heat
    double lnL

cdef struct Mc3SwapStats:
    uint64_t step
    uint64_t nswap

cdef enum:
    Mc3FreqProp             = 0
    Mc3RateProp             = 1
    Mc3RateShapeInvProp     = 2
    Mc3BrlenProp            = 3
    Mc3EtbrProp             = 4
    Mc3RateJumpProp         = 5
    Mc3PolytomyJumpProp     = 6
    Mc3RateShapeInvJumpProp = 7
    Mc3Prop                 = 8

cdef class Mc3Chain:
    cdef Mc3 master
    cdef unsigned run
    cdef unsigned ind
    cdef uint64_t nswap
    cdef double heat
    cdef unsigned swapInd
    cdef double swapProb
    cdef sfmt_t *swapPrng
    cdef sfmt_t *prng
    cdef Tree tree
    cdef Lik lik
    cdef double lnL
    cdef uint64_t step

    cdef bint freqPropose(self) except *
    cdef bint ratePropose(self) except *
    cdef bint rateShapeInvPropose(self) except *
    cdef bint brlenPropose(self) except *
    cdef bint etbrPropose(self) except *
    cdef void rateMergePropose(self, unsigned m0Ind, unsigned m1Ind, \
      list rclass, unsigned nrates) except *
    cdef void rateSplitPropose(self, unsigned m0Ind, unsigned m1Ind, \
      list rclass, unsigned nrates) except *
    cdef bint rateJumpPropose(self) except *
    cdef void polytomyMergePropose(self, Tree tree, unsigned nedges, \
      unsigned ntaxa) except *
    cdef void polytomySplitPropose(self, Tree tree, unsigned nedges, \
      unsigned ntaxa) except *
    cdef bint polytomyJumpPropose(self) except *
    cdef void rateShapeInvRemovePropose(self, unsigned m0Ind, unsigned m1Ind, \
      double alpha0) except *
    cdef void rateShapeInvAddPropose(self, unsigned m0Ind, unsigned m1Ind) \
      except *
    cdef bint rateShapeInvJumpPropose(self) except *
    cdef void advance(self) except *

cdef class Mc3:
    cdef readonly Alignment alignment
    cdef public str outPrefix
    cdef double _graphDelay

    # Convergence diagnostic parameters.
    cdef double _cvgDelay
    cdef double _cvgAlpha
    cdef double _cvgEpsilon

    # Chain parameters.
    cdef uint64_t _minStep
    cdef uint64_t _maxStep
    cdef unsigned _stride
    cdef unsigned _nruns

    # Metropolis coupling parameters.
    cdef unsigned _ncoupled
    cdef double _heatDelta
    cdef unsigned _swapStride

    # Gamma-distributed relative mutation rate parameters.
    cdef unsigned _ncat
    cdef bint _catMedian

    # Proposal parameters.
    cdef double _freqLambda
    cdef double _rateLambda
    cdef double _rateShapeInvLambda
    cdef double _brlenLambda
    cdef double _etbrPExt
    cdef double _etbrLambda

    # Model parameter priors.
    cdef double _rateShapeInvPrior
    cdef double _brlenPrior
    cdef double _rateJumpPrior
    cdef double _polytomyJumpPrior
    cdef double _rateShapeInvJumpPrior

    # Relative proposal probabilities, and corresponding CDF.
    cdef double props[8]
    cdef double propsCdf[9] # [1..8] corresponds to _*Prop fields.

    # Collated output files.
    cdef file tFile
    cdef file pFile
    cdef file sFile

    # Nested lists of chains.
    cdef list runs

    # Matrix of heat swapInfo structures, two for each pair of
    # Metropolis-coupled chains (even/odd steps).  The matrix is ordered as
    # such:
    #
    #                chain 0 chain 1 chain 2 chain 3
    #               |-------|-------|-------|-------|
    # run 0 chain 0 | ev,od |       |       |       |
    #       chain 1 |       |       |       |       |
    #       chain 2 |       |       |       |       |
    #       chain 3 |       |       |       |       |
    #               |-------|-------|-------|-------|
    # run 1 chain 0 |       |       |       |       |
    #       chain 1 |       |       |       |       |
    #       chain 2 |       |       |       |       |
    #       chain 3 |       |       |       |       |
    #               |-------|-------|-------|-------|
    cdef Mc3SwapInfo *swapInfo

    # Array of swap statistics, one element for each run.
    cdef Mc3SwapStats *swapStats

    # Nested lists of Lik samples from unheated chains.
#XXX    cdef list liks

    # Matrix of lnL samples from unheated chains.
    cdef double *lnLs
    cdef uint64_t lnLsMax

    # Scratch vector of sorted lnL samples, used by computeRcov().
    cdef double *rcovScratch
    cdef uint64_t rcovScratchMax

    cdef bint verbose

    cdef void sendSample(self, unsigned runInd, uint64_t step, double heat, \
      uint64_t nswap, Lik lik, double lnL) except *
    cdef void sendSwapInfo(self, unsigned runInd, unsigned srcChainInd, \
      unsigned dstChainInd, uint64_t step, double heat, double lnL) except *
    cdef void recvSwapInfo(self, unsigned runInd, unsigned dstChainInd, \
      unsigned srcChainInd, uint64_t step, double *heat, double *lnL) except *
    cdef double computeRcov(self, uint64_t last) except *
    cdef bint writeGraph(self, uint64_t sample, list rcovs) except *
    cdef str formatSwapStats(self, unsigned step)
    cpdef bint run(self, bint verbose=*) except *

    cdef double getGraphDelay(self)
    cdef void setGraphDelay(self, double graphDelay)
    # property graphDelay
    cdef double getCvgDelay(self)
    cdef void setCvgDelay(self, double cvgDelay)
    # property cvgDelay
    cdef double getCvgAlpha(self)
    cdef void setCvgAlpha(self, double cvgAlpha) except *
    # property cvgAlpha
    cdef double getCvgEpsilon(self)
    cdef void setCvgEpsilon(self, double cvgEpsilon) except *
    # property cvgEpsilon
    cdef uint64_t getMinStep(self)
    cdef void setMinStep(self, uint64_t minStep) except *
    # property minStep
    cdef uint64_t getMaxStep(self)
    cdef void setMaxStep(self, uint64_t maxStep) except *
    # property maxStep
    cdef unsigned getStride(self)
    cdef void setStride(self, unsigned stride) except *
    # property stride
    cdef unsigned getNruns(self)
    cdef void setNruns(self, unsigned nruns) except *
    # property nruns
    cdef unsigned getNcoupled(self)
    cdef void setNcoupled(self, unsigned ncoupled) except *
    # property ncoupled
    cdef double getHeatDelta(self)
    cdef void setHeatDelta(self, double heatDelta) except *
    # property heatDelta
    cdef unsigned getSwapStride(self)
    cdef void setSwapStride(self, unsigned swapStride) except *
    # property swapStride
    cdef unsigned getNcat(self)
    cdef void setNcat(self, unsigned ncat) except *
    # property ncat
    cdef bint getCatMedian(self)
    cdef void setCatMedian(self, bint catMedian)
    # property catMedian
    cdef double getFreqLambda(self)
    cdef void setFreqLambda(self, double freqLambda) except *
    # property freqLambda
    cdef double getRateLambda(self)
    cdef void setRateLambda(self, double rateLambda) except *
    # property rateLambda
    cdef double getRateShapeInvLambda(self)
    cdef void setRateShapeInvLambda(self, double rateShapeInvLambda) except *
    # property rateShapeInvLambda
    cdef double getBrlenLambda(self)
    cdef void setBrlenLambda(self, double brlenLambda) except *
    # property brlenLambda
    cdef double getEtbrPExt(self)
    cdef void setEtbrPExt(self, double etbrPExt) except *
    # property etbrPExt
    cdef double getEtbrLambda(self)
    cdef void setEtbrLambda(self, double etbrLambda) except *
    # property etbrLambda
    cdef double getRateShapeInvPrior(self)
    cdef void setRateShapeInvPrior(self, double rateShapeInvPrior) except *
    # property rateShapeInvPrior
    cdef double getBrlenPrior(self)
    cdef void setBrlenPrior(self, double brlenPrior) except *
    # property brlenPrior
    cdef double getRateJumpPrior(self)
    cdef void setRateJumpPrior(self, double rateJumpPrior) except *
    # property rateJumpPrior
    cdef double getPolytomyJumpPrior(self)
    cdef void setPolytomyJumpPrior(self, double polytomyJumpPrior) except *
    # property polytomyJumpPrior
    cdef double getRateShapeInvJumpPrior(self)
    cdef void setRateShapeInvJumpPrior(self, double rateShapeInvJumpPrior) \
      except *
    # property rateShapeInvJumpPrior
    cdef double getFreqProp(self)
    cdef void setFreqProp(self, double freqProp) except *
    # property freqProp
    cdef double getRateProp(self)
    cdef void setRateProp(self, double rateProp) except *
    # property rateProp
    cdef double getRateShapeInvProp(self)
    cdef void setRateShapeInvProp(self, double rateShapeInvProp) except *
    # property rateShapeInvProp
    cdef double getBrlenProp(self)
    cdef void setBrlenProp(self, double brlenProp) except *
    # property brlenProp
    cdef double getEtbrProp(self)
    cdef void setEtbrProp(self, double etbrProp) except *
    # property etbrProp
    cdef double getRateJumpProp(self)
    cdef void setRateJumpProp(self, double rateJumpProp) except *
    # property rateJumpProp
    cdef double getPolytomyJumpProp(self)
    cdef void setPolytomyJumpProp(self, double polytomyJumpProp) except *
    # property polytomyJumpProp
    cdef double getRateShapeInvJumpProp(self)
    cdef void setRateShapeInvJumpProp(self, double rateShapeInvJumpProp) \
      except *
    # property rateShapeInvJumpProp
