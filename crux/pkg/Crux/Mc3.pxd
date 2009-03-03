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

cdef struct Mc3RateStats:
    double ema # Exponential moving average.
    uint64_t n # Numerator.
    uint64_t d # Denominator.

cdef enum:
    Mc3WeightProp           =  0
    Mc3FreqProp             =  1
    Mc3RmultProp            =  2
    Mc3RateProp             =  3
    Mc3RateShapeInvProp     =  4
    Mc3BrlenProp            =  5
    Mc3EtbrProp             =  6
    Mc3RateJumpProp         =  7
    Mc3PolytomyJumpProp     =  8
    Mc3RateShapeInvJumpProp =  9
    Mc3Prop                 = 10

cdef class Mc3Chain:
    cdef Mc3 master
    cdef unsigned run
    cdef unsigned ind
    cdef uint64_t nswap
    cdef uint64_t accepts[Mc3Prop]
    cdef uint64_t rejects[Mc3Prop]
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
    cdef bint freqPropose(self) except *
    cdef bint rmultPropose(self) except *
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

    # Control parameters for diagnostic output.
    cdef double _graphDelay
    cdef double _emaAlpha

    # Convergence diagnostic parameters.
    cdef uint64_t _cvgSampStride
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

    # Mixture model parameters.
    cdef unsigned _nmodels

    # Gamma-distributed relative mutation rate parameters.
    cdef unsigned _ncat
    cdef bint _catMedian

    # Proposal parameters.
    cdef double _weightLambda
    cdef double _freqLambda
    cdef double _rmultLambda
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
    cdef double props[Mc3Prop]
    cdef double propsCdf[Mc3Prop+1] # [1..Mc3Prop] corresponds to props.

    # Output files.
    cdef file lFile
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
    cdef Mc3RateStats *swapStats

    # Arrays of accept/reject statistics, one element for each run.
    cdef Mc3RateStats *propStats[Mc3Prop]

    # Array of pointers to lnL sample arrays from unheated chains, plus one
    # extra (used as a scratch area).
    cdef double **lnLs
    cdef uint64_t lnLsMax

    cdef bint verbose

    cdef void sendSample(self, unsigned runInd, uint64_t step, double heat, \
      uint64_t nswap, uint64_t *accepts, uint64_t *rejects, Lik lik, \
      double lnL) except *
    cdef void sendSwapInfo(self, unsigned runInd, unsigned srcChainInd, \
      unsigned dstChainInd, uint64_t step, double heat, double lnL) except *
    cdef void recvSwapInfo(self, unsigned runInd, unsigned dstChainInd, \
      unsigned srcChainInd, uint64_t step, double *heat, double *lnL) except *
    cdef void initLogs(self) except *
    cdef double computeRcov(self, uint64_t last) except *
    cdef bint writeGraph(self, uint64_t sample) except *
    cdef str formatLnLs(self, uint64_t sample, str fmt)
    cdef str formatRateStats(self, Mc3RateStats *rateStats)
    cdef str formatPropStats(self)
    cdef void updateDiags(self, uint64_t step)
    cpdef bint run(self, bint verbose=*) except *

    cdef double getGraphDelay(self)
    cdef void setGraphDelay(self, double graphDelay)
    # property graphDelay
    cdef double getEmaAlpha(self)
    cdef void setEmaAlpha(self, double emaAlpha) except *
    # property emaAlpha
    cdef uint64_t getCvgSampStride(self)
    cdef void setCvgSampStride(self, uint64_t cvgSampStride)
    # property cvgSampStride
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
    cdef unsigned getNmodels(self)
    cdef void setNmodels(self, unsigned nmodels) except *
    # property nmodels
    cdef unsigned getNcat(self)
    cdef void setNcat(self, unsigned ncat) except *
    # property ncat
    cdef bint getCatMedian(self)
    cdef void setCatMedian(self, bint catMedian)
    # property catMedian
    cdef double getWeightLambda(self)
    cdef void setWeightLambda(self, double weightLambda) except *
    # property weightLambda
    cdef double getFreqLambda(self)
    cdef void setFreqLambda(self, double freqLambda) except *
    # property freqLambda
    cdef double getRmultLambda(self)
    cdef void setRmultLambda(self, double rmultLambda) except *
    # property rmultLambda
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
    cdef double getWeightProp(self)
    cdef void setWeightProp(self, double weightProp) except *
    # property weightProp
    cdef double getFreqProp(self)
    cdef void setFreqProp(self, double freqProp) except *
    # property freqProp
    cdef double getRmultProp(self)
    cdef void setRmultProp(self, double rmultProp) except *
    # property rmultProp
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
