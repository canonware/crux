/******************************************************************************
 *
 * <Copyright = jasone>
 * <License>
 *
 ******************************************************************************
 *
 * Version: Crux <Version = crux>
 *
 ******************************************************************************/

void *
mem_malloc_e(size_t a_size, const char *a_filename, uint32_t a_line_num);

void *
mem_calloc_e(size_t a_number, size_t a_size, const char *a_filename,
	     uint32_t a_line_num);

void *
mem_realloc_e(void *a_ptr, size_t a_size, const char *a_filename,
	      uint32_t a_line_num);

void
mem_free_e(void *a_ptr, const char *a_filename, uint32_t a_line_num);

/* These macros are declared differently, depending on whether this is a debug
 * library, because consistently using arguments of NULL and 0 reduces the size
 * of the generated binary.  Since these arguments aren't used in the non-debug
 * version anyway, this is a free (though perhaps small) memory savings. */
#ifdef CW_DBG
#define cw_malloc(a_size)						\
    mem_malloc_e((a_size), __FILE__, __LINE__)
#define cw_calloc(a_number, a_size)					\
    mem_calloc_e((a_number), (a_size), __FILE__, __LINE__)
#define cw_realloc(a_ptr, a_size)					\
    mem_realloc_e((a_ptr), (a_size), __FILE__, __LINE__)
#define cw_free(a_ptr)							\
    mem_free_e((a_ptr), __FILE__, __LINE__)

#else
#define cw_malloc(a_size)						\
    mem_malloc_e((a_size), NULL, 0)
#define cw_calloc(a_number, a_size)					\
    mem_calloc_e((a_number), (a_size), NULL, 0)
#define cw_realloc(a_ptr, a_size)					\
    mem_realloc_e((a_ptr), (a_size), NULL, 0)
#define cw_free(a_ptr)							\
    mem_free_e((a_ptr), NULL, 0)

#endif
