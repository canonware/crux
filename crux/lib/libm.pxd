cdef extern from "math.h":
    cdef double HUGE_VALF

    cdef double log(double x)

cdef extern from "float.h":
    cdef float FLT_EPSILON
