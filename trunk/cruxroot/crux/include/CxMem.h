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
CxMemMallocE(size_t aSize, const char *aFilename, uint32_t aLineNum);

void *
CxMemCallocE(size_t aNumber, size_t aSize, const char *aFilename,
	     uint32_t aLineNum);

void *
CxMemReallocE(void *aPtr, size_t aSize, const char *aFilename,
	      uint32_t aLineNum);

void
CxMemFreeE(void *aPtr, const char *aFilename, uint32_t aLineNum);

/* These macros are declared differently, depending on whether this is a debug
 * library, because consistently using arguments of NULL and 0 reduces the size
 * of the generated binary.  Since these arguments aren't used in the non-debug
 * version anyway, this is a free (though perhaps small) memory savings. */
#ifdef CxmDebug
#define CxmMalloc(aSize)						\
    CxMemMallocE((aSize), __FILE__, __LINE__)
#define CxmCalloc(aNumber, aSize)					\
    CxMemCallocE((aNumber), (aSize), __FILE__, __LINE__)
#define CxmRealloc(aPtr, aSize)						\
    CxMemReallocE((aPtr), (aSize), __FILE__, __LINE__)
#define CxmFree(aPtr)							\
    CxMemFreeE((aPtr), __FILE__, __LINE__)

#else
#define CxmMalloc(aSize)						\
    CxMemMallocE((aSize), NULL, 0)
#define CxmCalloc(aNumber, aSize)					\
    CxMemCallocE((aNumber), (aSize), NULL, 0)
#define CxmRealloc(aPtr, aSize)						\
    CxMemReallocE((aPtr), (aSize), NULL, 0)
#define CxmFree(aPtr)							\
    CxMemFreeE((aPtr), NULL, 0)

#endif
