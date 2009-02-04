from libc cimport uint32_t, uint64_t

cdef extern from "SFMT.h":
    ctypedef void sfmt_t

    cdef uint32_t gen_rand32(sfmt_t *ctx)
    cdef uint32_t gen_rand32_range(sfmt_t *ctx, uint32_t limit)
    cdef uint64_t gen_rand64(sfmt_t *ctx)
    cdef uint64_t gen_rand64_range(sfmt_t *ctx, uint64_t limit)
    cdef void fill_array32(sfmt_t *ctx, uint32_t *array, int size)
    cdef void fill_array64(sfmt_t *ctx, uint64_t *array, int size)
    cdef sfmt_t *init_gen_rand(uint32_t seed)
    cdef sfmt_t *init_by_array(uint32_t *init_key, int key_length)
    cdef void fini_gen_rand(sfmt_t *ctx)
    cdef char *get_idstring()
    cdef int get_min_array_size32()
    cdef int get_min_array_size64()
    cdef inline double to_real1(uint32_t v)
    cdef inline double genrand_real1(sfmt_t *ctx)
    cdef inline double to_real2(uint32_t v)
    cdef inline double genrand_real2(sfmt_t *ctx)
    cdef inline double to_real3(uint32_t v)
    cdef inline double genrand_real3(sfmt_t *ctx)
    cdef inline double to_res53(uint64_t v)
    cdef inline double to_res53_mix(uint32_t x, uint32_t y)
    cdef inline double genrand_res53(sfmt_t *ctx)
    cdef inline double genrand_res53_mix(sfmt_t *ctx)
