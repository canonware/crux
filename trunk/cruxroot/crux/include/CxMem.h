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
CxMemMallocE(size_t a_size, const char *a_filename, uint32_t a_line_num);

void *
CxMemCallocE(size_t a_number, size_t a_size, const char *a_filename,
	     uint32_t a_line_num);

void *
CxMemReallocE(void *a_ptr, size_t a_size, const char *a_filename,
	      uint32_t a_line_num);

void
CxMemFreeE(void *a_ptr, const char *a_filename, uint32_t a_line_num);

/* These macros are declared differently, depending on whether this is a debug
 * library, because consistently using arguments of NULL and 0 reduces the size
 * of the generated binary.  Since these arguments aren't used in the non-debug
 * version anyway, this is a free (though perhaps small) memory savings. */
#ifdef CxmDebug
#define CxmMalloc(a_size)						\
    CxMemMallocE((a_size), __FILE__, __LINE__)
#define CxmCalloc(a_number, a_size)					\
    CxMemCallocE((a_number), (a_size), __FILE__, __LINE__)
#define CxmRealloc(a_ptr, a_size)					\
    CxMemReallocE((a_ptr), (a_size), __FILE__, __LINE__)
#define CxmFree(a_ptr)							\
    CxMemFreeE((a_ptr), __FILE__, __LINE__)

#else
#define CxmMalloc(a_size)						\
    CxMemMallocE((a_size), NULL, 0)
#define CxmCalloc(a_number, a_size)					\
    CxMemCallocE((a_number), (a_size), NULL, 0)
#define CxmRealloc(a_ptr, a_size)					\
    CxMemReallocE((a_ptr), (a_size), NULL, 0)
#define CxmFree(a_ptr)							\
    CxMemFreeE((a_ptr), NULL, 0)

#endif
