cdef extern from "math.h":
    cdef double HUGE_VALF
    cdef double NAN
    cdef double INFINITY

    cdef double exp(double x)
    cdef double log(double x)
    cdef float logf(float x)
    cdef double pow(double x, double y)
    cdef double fabs(double x)
    cdef float roundf(float x)
    cdef double sqrt(double x)

cdef extern from "float.h":
    cdef float FLT_EPSILON
