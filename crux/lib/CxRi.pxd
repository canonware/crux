from libc cimport uint32_t
from SFMT cimport sfmt_t

cdef extern from "CxRi.h":
    ctypedef struct CxtRi:
        sfmt_t *prng
        uint32_t arr
        uint32_t arrLen
        uint32_t nints
        uint32_t ind

    cdef void CxRiNew(CxtRi *aRi, sfmt_t *aPrng)
    cdef void CxRiDelete(CxtRi *aRi)
    cdef bint CxRiInit(CxtRi *aRi, uint32_t aNints)
    cdef uint32_t CxRiNintsGet(CxtRi *aRi)
    cdef uint32_t CxRiIndGet(CxtRi *aRi)
    cdef uint32_t CxRiRandomGet(CxtRi *aRi)
