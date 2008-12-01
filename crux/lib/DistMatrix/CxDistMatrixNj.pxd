from CxDistMatrix cimport *

cdef extern from "CxDistMatrixNj.h":
    cdef inline int CxDistMatrixNjDistCompare(CxtDMDist a, CxtDMDist b)
