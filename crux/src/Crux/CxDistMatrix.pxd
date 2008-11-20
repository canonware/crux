cdef extern from "CxDistMatrix.h":
    ctypedef double CxtDMDist
    ctypedef unsigned long CxtDMSize

    cdef inline CxtDMSize CxDistMatrixNxy2i(CxtDMSize n, CxtDMSize x,
      CxtDMSize y)
    cdef CxtDMDist *CxDistMatrixNew(CxtDMSize ntaxa)
    cdef void CxDistMatrixDelete(CxtDMDist *matrix)
