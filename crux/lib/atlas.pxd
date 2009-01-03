cdef extern from "cblas.h":
    cdef enum CBLAS_ORDER:
        CblasRowMajor

    cdef void cblas_dscal(int N, double alpha, double *X, int incX)
    cdef double cblas_dasum(int N, double *X, int incX)

cdef extern from "clapack.h":
    cdef int clapack_dgetrf(CBLAS_ORDER Order, int M, int N, double *A, int lda,
      int *ipiv)
