cdef extern from "stdint.h":
    ctypedef unsigned uint32_t

cdef extern from "CxRi.h":
    ctypedef struct CxtRi:
        uint32_t arr
        uint32_t arrLen
        uint32_t nints
        uint32_t ind

    cdef void CxRiNew(CxtRi *aRi)
    cdef void CxRiDelete(CxtRi *aRi)
    cdef bint CxRiInit(CxtRi *aRi, uint32_t aNints)
    cdef uint32_t CxRiNintsGet(CxtRi *aRi)
    cdef uint32_t CxRiIndGet(CxtRi *aRi)
    cdef uint32_t CxRiRandomGet(CxtRi *aRi)
