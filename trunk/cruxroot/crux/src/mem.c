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
mem_malloc_e(size_t a_size, const char *a_filename, uint32_t a_line_num)
{
    void *retval;

    retval = malloc(a_size);
    if (retval == NULL)
    {
#ifdef CW_DBG
	fprintf(stderr, "%s(): %p <-- malloc(%zu) at %s:%u\n",
		__func__, retval, a_size, a_filename, a_line_num);
#endif
	xep_throw(CW_CRUXX_OOM);
    }

    return retval;
}

void *
mem_calloc_e(size_t a_number, size_t a_size, const char *a_filename,
	     uint32_t a_line_num)
{
    void *retval;

    retval = calloc(a_number, a_size);
    if (retval == NULL)
    {
#ifdef CW_DBG
	fprintf(stderr, "%s(): %p <-- calloc(%zu, %zu) at %s:%u\n",
		__func__, retval, a_number, a_size,
		a_filename, a_line_num);
#endif
	xep_throw(CW_CRUXX_OOM);
    }
    return retval;
}

void *
mem_realloc_e(void *a_ptr, size_t a_size, const char *a_filename,
	      uint32_t a_line_num)
{
    void *retval;

    retval = realloc(a_ptr, a_size);
    if (retval == NULL)
    {
#ifdef CW_DBG
	fprintf(stderr, "%s(): %p <-- realloc(%p, %zu) at %s:%u\n",
		__func__, retval, a_ptr, a_size,
		a_filename, a_line_num);
#endif
	xep_throw(CW_CRUXX_OOM);
    }

    return retval;
}

void
mem_free_e(void *a_ptr, const char *a_filename, uint32_t a_line_num)
{
    free(a_ptr);
}
