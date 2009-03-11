from libc cimport *
from Crux.Mc3 cimport Mc3
from Crux.Tree cimport Tree

cdef class Msamp:
    cdef double weight
    cdef double rmult
    cdef str rclass
    cdef list rates
    cdef double wNorm
    cdef double alpha
    cdef list freqs

cdef class Samp:
    cdef double lnL  # Call Post.parseS() before accessing.
    cdef list msamps # Call Post.parseP() before accessing.
    cdef Tree tree   # Call Post.parseT() before accessing.

cdef class Post:
    cdef bint verbose
    cdef readonly Mc3 mc3
    cdef readonly uint64_t seed
    cdef readonly uint64_t burnin
    cdef readonly uint64_t nsamples
    cdef readonly uint64_t stepFirst
    cdef readonly uint64_t stepLast
    cdef bint _sDone, _pDone, _tDone
    cdef readonly unsigned maxModels # Call parseP() before accessing.

    cdef double _rcovLnL
    cdef double _rcovRclass

    cdef double _lnLMinCred, _lnLMaxCred, _lnLMean, _lnLVar, _lnLMed
    cdef double _nmodelsMinCred, _nmodelsMaxCred, _nmodelsMean
    cdef double _nmodelsVar, _nmodelsMed

    cdef list runs

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
    # property lnLMinCred
    cdef double getNmodelsMaxCred(self) except -1.0
    # property lnLMaxCred
    cdef double getNmodelsMean(self) except -1.0
    # property lnLMean
    cdef double getNmodelsVar(self) except -1.0
    # property lnLVar
    cdef double getNmodelsMed(self) except -1.0
    # property lnLMed
