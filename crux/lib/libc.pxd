cdef extern from "sys/types.h":
    ctypedef unsigned long size_t

cdef extern from "stdint.h":
    ctypedef int int32_t
    ctypedef unsigned uint32_t

cdef extern from "stdlib.h":
    cdef void *malloc(size_t size)
    cdef void *calloc(size_t nmemb, size_t size)
    cdef void free(void *ptr)
    cdef int abs(int i)

cdef extern from "string.h":
    cdef void *memset(void *s, int c, size_t n)
    cdef void *memcpy(void *dest, void *src, size_t n)

cdef extern from "strings.h":
    cdef int ffs(int i)