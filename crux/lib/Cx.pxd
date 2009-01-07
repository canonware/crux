cdef extern from "Cx.h":
    cdef void CxInit()
    cdef inline int CxCmp2Richcmp(int cmp, int op)
