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
CxMemMallocE(size_t a_size, const char *a_filename, uint32_t a_line_num)
{
    void *retval;

    retval = malloc(a_size);
    if (retval == NULL)
    {
#ifdef CxmDebug
	fprintf(stderr, "%s(): %p <-- malloc(%zu) at %s:%u\n",
		__func__, retval, a_size, a_filename, a_line_num);
#endif
	CxmXepThrow(CxmXepOOM);
    }

    return retval;
}

void *
CxMemCallocE(size_t a_number, size_t a_size, const char *a_filename,
	     uint32_t a_line_num)
{
    void *retval;

    retval = calloc(a_number, a_size);
    if (retval == NULL)
    {
#ifdef CxmDebug
	fprintf(stderr, "%s(): %p <-- calloc(%zu, %zu) at %s:%u\n",
		__func__, retval, a_number, a_size,
		a_filename, a_line_num);
#endif
	CxmXepThrow(CxmXepOOM);
    }
    return retval;
}

void *
CxMemReallocE(void *a_ptr, size_t a_size, const char *a_filename,
	      uint32_t a_line_num)
{
    void *retval;

    retval = realloc(a_ptr, a_size);
    if (retval == NULL)
    {
#ifdef CxmDebug
	fprintf(stderr, "%s(): %p <-- realloc(%p, %zu) at %s:%u\n",
		__func__, retval, a_ptr, a_size,
		a_filename, a_line_num);
#endif
	CxmXepThrow(CxmXepOOM);
    }

    return retval;
}

void
CxMemFreeE(void *a_ptr, const char *a_filename, uint32_t a_line_num)
{
    free(a_ptr);
}
