cdef extern from "math.h":
    cdef double HUGE_VALF

    cdef double log(double x)
    cdef float logf(float x)
    cdef double pow(double x, double y)

    cdef float roundf(float x)

cdef extern from "float.h":
    cdef float FLT_EPSILON
