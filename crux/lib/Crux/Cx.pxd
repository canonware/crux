cdef extern from "Cx.h":
    cdef unsigned CxNcpus

    cdef void CxInit()
    cdef void CxThreaded()
    cdef inline int CxCmp2Richcmp(int cmp, int op)
