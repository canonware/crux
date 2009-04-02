from libc cimport *
from Crux.Mc3 cimport Mc3
from Crux.Tree cimport Tree
from Crux.Tree.Sumt cimport Sumt

cdef class Msamp:
    cdef readonly double weight
    cdef readonly double rmult
    cdef readonly str rclass
    cdef readonly list rates
    cdef readonly double wNorm
    cdef readonly double alpha
    cdef readonly list freqs

cdef class Samp:
    cdef readonly double lnL  # Call Post.parseS() before accessing.
    cdef readonly list msamps # Call Post.parseP() before accessing.
    cdef readonly Tree tree   # Call Post.parseT() before accessing.

cdef class Post:
    cdef bint verbose
    cdef readonly Mc3 mc3
    cdef readonly uint64_t seed
    cdef readonly uint64_t burnin
    cdef bint _sDone, _pDone, _tDone
    cdef readonly uint64_t nsamples  # Call parseS() before accessing.
    cdef readonly uint64_t stepFirst # Call parseS() before accessing.
    cdef readonly uint64_t stepLast  # Call parseS() before accessing.
    cdef readonly unsigned maxModels # Call parseP() before accessing.

    cdef double _rcovLnL
    cdef double _rcovRclass

    cdef double _lnLMinCred, _lnLMaxCred, _lnLMean, _lnLVar, _lnLMed
    cdef double _nmodelsMinCred, _nmodelsMaxCred, _nmodelsMean
    cdef double _nmodelsVar, _nmodelsMed
    cdef list _ratesMinCred, _ratesMaxCred, _ratesMean, _ratesVar, _ratesMed
    cdef list _freqsMinCred, _freqsMaxCred, _freqsMean, _freqsVar, _freqsMed
    cdef double _alphaMinCred, _alphaMaxCred, _alphaMean, _alphaVar, _alphaMed
    cdef Sumt _sumt

    cdef readonly list runs

    cdef void _parseL(self) except *
    cpdef parseS(self)
    cpdef parseP(self)
    cpdef parseT(self)

    cdef void _computeRcovLnL(self) except *
    cdef double getRcovLnL(self) except -1.0
    # property rcovLnL
    cdef void _computeRcovRclass(self) except *
    cdef double getRcovRclass(self) except -1.0
    # property rcovRclass

    cdef void _summarizeLnL(self) except *
    cdef double getLnLMinCred(self) except -1.0
    # property lnLMinCred
    cdef double getLnLMaxCred(self) except -1.0
    # property lnLMaxCred
    cdef double getLnLMean(self) except -1.0
    # property lnLMean
    cdef double getLnLVar(self) except -1.0
    # property lnLVar
    cdef double getLnLMed(self) except -1.0
    # property lnLMed

    cdef void _summarizeNmodels(self) except *
    cdef double getNmodelsMinCred(self) except -1.0
    # property nmodelsMinCred
    cdef double getNmodelsMaxCred(self) except -1.0
    # property nmodelsMaxCred
    cdef double getNmodelsMean(self) except -1.0
    # property nmodelsMean
    cdef double getNmodelsVar(self) except -1.0
    # property nmodelsVar
    cdef double getNmodelsMed(self) except -1.0
    # property nmodelsMed

    cdef void _summarizeRates(self) except *
    cdef list getRatesMinCred(self)
    # property ratesMinCred
    cdef list getRatesMaxCred(self)
    # property ratesMaxCred
    cdef list getRatesMean(self)
    # property ratesMean
    cdef list getRatesVar(self)
    # property ratesVar
    cdef list getRatesMed(self)
    # property ratesMed

    cdef void _summarizeFreqs(self) except *
    cdef list getFreqsMinCred(self)
    # property freqsMinCred
    cdef list getFreqsMaxCred(self)
    # property freqsMaxCred
    cdef list getFreqsMean(self)
    # property freqsMean
    cdef list getFreqsVar(self)
    # property freqsVar
    cdef list getFreqsMed(self)
    # property freqsMed

    cdef void _summarizeAlpha(self) except *
    cdef double getAlphaMinCred(self) except -1.0
    # property alphaMinCred
    cdef double getAlphaMaxCred(self) except -1.0
    # property alphaMaxCred
    cdef double getAlphaMean(self) except -1.0
    # property alphaMean
    cdef double getAlphaVar(self) except -1.0
    # property alphaVar
    cdef double getAlphaMed(self) except -1.0
    # property alphaMed

    cdef void _summarizeTrees(self) except *
    cdef Sumt getSumt(self)
    # property sumt
