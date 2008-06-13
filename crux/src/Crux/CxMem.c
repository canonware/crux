//==============================================================================
//
// <Copyright = jasone>
// <License>
//
//==============================================================================
//
// Version: Crux <Version = crux>
//
//==============================================================================

#include "Crux/_cruxmodule.h"

void *
CxMemMallocE(size_t aSize, const char *aFilename, uint32_t aLineNum)
{
    void *rVal;

    rVal = malloc(aSize);
    if (rVal == NULL)
    {
#ifdef CxmDebug
	fprintf(stderr, "%s(): %p <-- malloc(%zu) at %s:%u\n",
		__func__, rVal, aSize, aFilename, aLineNum);
#endif
	CxmXepThrow(CxmXepOOM);
    }

    return rVal;
}

void *
CxMemCallocE(size_t aNumber, size_t aSize, const char *aFilename,
	     uint32_t aLineNum)
{
    void *rVal;

    rVal = calloc(aNumber, aSize);
    if (rVal == NULL)
    {
#ifdef CxmDebug
	fprintf(stderr, "%s(): %p <-- calloc(%zu, %zu) at %s:%u\n",
		__func__, rVal, aNumber, aSize,
		aFilename, aLineNum);
#endif
	CxmXepThrow(CxmXepOOM);
    }
    return rVal;
}

void *
CxMemReallocE(void *aPtr, size_t aSize, const char *aFilename,
	      uint32_t aLineNum)
{
    void *rVal;

    rVal = realloc(aPtr, aSize);
    if (rVal == NULL)
    {
#ifdef CxmDebug
	fprintf(stderr, "%s(): %p <-- realloc(%p, %zu) at %s:%u\n",
		__func__, rVal, aPtr, aSize,
		aFilename, aLineNum);
#endif
	CxmXepThrow(CxmXepOOM);
    }

    return rVal;
}

void
CxMemFreeE(void *aPtr, const char *aFilename, uint32_t aLineNum)
{
    free(aPtr);
}
