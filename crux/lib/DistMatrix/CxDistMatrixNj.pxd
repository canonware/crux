from CxDistMatrix cimport *

cdef extern from "CxDistMatrixNj.h":
    cdef inline int CxDistMatrixNjDistCompare(CxtDMDist a, CxtDMDist b)
    cdef inline int CxDistMatrixNjDistLt(CxtDMDist a, CxtDMDist b)
    cdef inline int CxDistMatrixNjDistEq(CxtDMDist a, CxtDMDist b)
