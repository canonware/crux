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

#include "../include/_cruxmodule.h"

void *
CxMemMallocE(size_t aSize, const char *aFilename, uint32_t aLineNum)
{
    void *retval;

    retval = malloc(aSize);
    if (retval == NULL)
    {
#ifdef CxmDebug
	fprintf(stderr, "%s(): %p <-- malloc(%zu) at %s:%u\n",
		__func__, retval, aSize, aFilename, aLineNum);
#endif
	CxmXepThrow(CxmXepOOM);
    }

    return retval;
}

void *
CxMemCallocE(size_t aNumber, size_t aSize, const char *aFilename,
	     uint32_t aLineNum)
{
    void *retval;

    retval = calloc(aNumber, aSize);
    if (retval == NULL)
    {
#ifdef CxmDebug
	fprintf(stderr, "%s(): %p <-- calloc(%zu, %zu) at %s:%u\n",
		__func__, retval, aNumber, aSize,
		aFilename, aLineNum);
#endif
	CxmXepThrow(CxmXepOOM);
    }
    return retval;
}

void *
CxMemReallocE(void *aPtr, size_t aSize, const char *aFilename,
	      uint32_t aLineNum)
{
    void *retval;

    retval = realloc(aPtr, aSize);
    if (retval == NULL)
    {
#ifdef CxmDebug
	fprintf(stderr, "%s(): %p <-- realloc(%p, %zu) at %s:%u\n",
		__func__, retval, aPtr, aSize,
		aFilename, aLineNum);
#endif
	CxmXepThrow(CxmXepOOM);
    }

    return retval;
}

void
CxMemFreeE(void *aPtr, const char *aFilename, uint32_t aLineNum)
{
    free(aPtr);
}
