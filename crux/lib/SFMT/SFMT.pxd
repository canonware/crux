cdef extern from "stdint.h":
    ctypedef unsigned uint32_t
    ctypedef unsigned long long uint64_t

cdef extern from "SFMT.h":
    cdef uint32_t gen_rand32()
    cdef uint32_t gen_rand32_range(uint32_t range)
    cdef uint64_t gen_rand64()
    cdef uint64_t gen_rand64_range(uint64_t range)
    cdef void fill_array32(uint32_t *array, int size)
    cdef void fill_array64(uint64_t *array, int size)
    cdef void init_gen_rand(uint32_t seed)
    cdef void init_by_array(uint32_t *init_key, int key_length)
    cdef char *get_idstring()
    cdef int get_min_array_size32()
    cdef int get_min_array_size64()
    cdef inline double to_real1(uint32_t v)
    cdef inline double genrand_real1()
    cdef inline double to_real2(uint32_t v)
    cdef inline double genrand_real2()
    cdef inline double to_real3(uint32_t v)
    cdef inline double genrand_real3()
    cdef inline double to_res53(uint64_t v)
    cdef inline double to_res53_mix(uint32_t x, uint32_t y)
    cdef inline double genrand_res53()
    cdef inline double genrand_res53_mix()
