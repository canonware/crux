# Forward declaration.
cdef class Mc3

from libc cimport uint64_t
from SFMT cimport *
from Crux.CTMatrix cimport Alignment
from Crux.Mc3.Chain cimport Chain, PropCnt
from Crux.Tree cimport Tree
from Crux.Tree.Lik cimport Lik
IF @enable_mpi@:
    cimport mpi4py.mpi_c as mpi
from Crux.Mc3.Post cimport Post

cdef struct Mc3SwapInfo:
    uint64_t step
    double heat
    double lnL

cdef struct Mc3RateStats:
    double ema # Exponential moving average.
    uint64_t n # Numerator.
    uint64_t d # Denominator.

cdef class Mc3:
    cdef readonly Alignment alignment
    cdef public str outPrefix

    # If set, used as the basis for computing topological risk.
    cdef Post _prelim
    cdef list _prelimParts

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

    # Use +I models if true.
    cdef bint _invar

    # Proposal parameters.
    cdef double _weightLambda
    cdef double _freqLambda
    cdef double _rmultLambda
    cdef double _rateLambda
    cdef double _rateShapeInvLambda
    cdef double _invarLambda
    cdef double _brlenLambda
    cdef double _etbrPExt
    cdef double _etbrLambda

    # Model parameter priors.
    cdef double _rateShapeInvPrior
    cdef double _invarPrior
    cdef double _brlenPrior
    cdef double _rateJumpPrior
    cdef double _polytomyJumpPrior
    cdef double _rateShapeInvJumpPrior
    cdef double _invarJumpPrior
    cdef double _freqJumpPrior
    cdef double _mixtureJumpPrior

    # Relative proposal probabilities, and corresponding CDF.
    cdef double props[PropCnt]
    cdef double propsPdf[PropCnt] # Sums to 1.0.
    cdef double propsCdf[PropCnt+1] # [1..PropCnt] corresponds to props.

    # Output files.
    cdef file lFile
    cdef file tFile
    cdef file pFile
    cdef file sFile

    # Nested lists of chains.
    cdef list runs

    IF @enable_mpi@:
        cdef int mpiPpn
        cdef int mpiSize
        cdef int mpiRank
        cdef bint mpiActive
        cdef bint mpiActiveCommAlloced
        cdef mpi.MPI_Comm mpiActiveComm

    # Matrix of heat swapInfo structures, two for each pair of
    # Metropolis-coupled chains (even/odd swaps).  The matrix is ordered as
    # such:
    #
    #             \to
    #          from\ chain 0 chain 1 chain 2 chain 3
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
    cdef Mc3RateStats *propStats[PropCnt]

    # Temporary storage for the most recent Liks.  These are used to write
    # statistics to pFile/tFile.
    cdef list liks

    # Temporary storage for the most recent lnLs.  These are appended to the
    # lnLs arrays once every run has provided its lnL via sendSample().
    cdef double *cachedLnLs

    # Array of pointers to lnL sample arrays from unheated chains, plus one
    # extra (used as a scratch area).
    cdef double **lnLs
    cdef size_t lnLsMax

    cdef bint verbose

    cpdef Mc3 dup(self)
    IF @enable_mpi@:
        cdef void storeSwapStats(self) except *
        cdef void storePropStats(self) except *
    cdef void updateDiags(self, uint64_t step) except *
    cdef void storeLiksLnLsUni(self, uint64_t step) except *
    IF @enable_mpi@:
        cdef void storeLiksLnLsMpi(self, uint64_t step) except *
    cdef void storeLiksLnLs(self, uint64_t step) except *
    cdef void sendSample(self, unsigned runInd, uint64_t step, double heat, \
      uint64_t nswap, uint64_t *accepts, uint64_t *rejects, Lik lik, \
      double lnL) except *
    cdef void sendSwapInfo(self, unsigned runInd, unsigned srcChainInd, \
      unsigned dstChainInd, uint64_t step, double heat, double lnL) except *
    cdef void recvSwapInfoUni(self, unsigned runInd, unsigned dstChainInd, \
      unsigned srcChainInd, uint64_t step, double *heat, double *lnL) except *
    IF @enable_mpi@:
        cdef void recvSwapInfoMpi(self, unsigned runInd, unsigned dstChainInd, \
          unsigned srcChainInd, uint64_t step, double *heat, double *lnL) \
          except *
    cdef void recvSwapInfo(self, unsigned runInd, unsigned dstChainInd, \
      unsigned srcChainInd, uint64_t step, double *heat, double *lnL) except *
    IF @enable_mpi@:
        cdef list rank2interleave(self, unsigned nnodes)
        cdef bint initActive(self) except *
    cdef void initLogs(self) except *
    cdef void initPrelim(self) except *
    cdef void initSwapInfo(self) except *
    cdef void initSwapStats(self) except *
    cdef void initPropStats(self) except *
    cdef void initLiks(self) except *
    cdef void resetLiks(self) except *
    cdef void initLnLs(self) except *
    cdef void initProps(self) except *
    cdef void initRunsUni(self, list liks) except *
    IF @enable_mpi@:
        cdef void initRunsMpi(self, list liks) except *
    cdef void initRuns(self, list liks) except *
    cdef void advanceUni(self) except *
    IF @enable_mpi@:
        cdef void advanceMpi(self) except *
    cdef void advance(self) except *
    cdef double computeRcovUni(self, uint64_t step) except *
    IF @enable_mpi@:
        cdef double computeRcovMpi(self, uint64_t step) except *
    cdef double computeRcov(self, uint64_t step) except *
    cdef bint writeGraph(self, uint64_t step) except *
    cdef str formatRclass(self, Lik lik, unsigned model)
    cdef str formatRates(self, Lik lik, unsigned model, str fmt)
    cdef str formatFreqs(self, Lik lik, unsigned model, str fmt)
    cdef str formatLnLs(self, uint64_t step, str fmt)
    cdef str formatRateStats(self, Mc3RateStats *rateStats)
    cdef str formatPropStats(self)
    cdef void lWrite(self, str s) except *
    cdef void tWrite(self, uint64_t step) except *
    cdef void pWrite(self, uint64_t step) except *
    cdef void sWrite(self, uint64_t step, double rcov) except *
    cdef bint sample(self, uint64_t step) except *
    cdef void randomDnaQ(self, Lik lik, unsigned model, sfmt_t *prng) except *
    cpdef Lik randomLik(self, Tree tree=*)
    cpdef run(self, bint verbose=*, list liks=*)

    cdef Post getPrelim(self)
    cdef void setPrelim(self, Post prelim) except *
    # property prelim

    IF @enable_mpi@:
        cdef int getPpn(self)
        cdef void setPpn(self, int ppn) except *
        # property ppn

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
    cdef bint getInvar(self)
    cdef void setInvar(self, bint invar)
    # property invar
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
    cdef double getInvarLambda(self)
    cdef void setInvarLambda(self, double invarLambda) except *
    # property invarLambda
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
    cdef double getInvarPrior(self)
    cdef void setInvarPrior(self, double invarPrior) except *
    # property invarPrior
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
    cdef double getInvarJumpPrior(self)
    cdef void setInvarJumpPrior(self, double invarJumpPrior) except *
    # property invarJumpPrior
    cdef double getFreqJumpPrior(self)
    cdef void setFreqJumpPrior(self, double freqJumpPrior) except *
    # property freqJumpPrior
    cdef double getMixtureJumpPrior(self)
    cdef void setMixtureJumpPrior(self, double mixtureJumpPrior) except *
    # property mixtureJumpPrior
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
    cdef double getInvarProp(self)
    cdef void setInvarProp(self, double invarProp) except *
    # property invarProp
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
    cdef double getInvarJumpProp(self)
    cdef void setInvarJumpProp(self, double invarJumpProp) except *
    # property invarJumpProp
    cdef double getFreqJumpProp(self)
    cdef void setFreqJumpProp(self, double freqJumpProp) except *
    # property freqJumpProp
    cdef double getMixtureJumpProp(self)
    cdef void setMixtureJumpProp(self, double mixtureJumpProp) except *
    # property mixtureJumpProp
