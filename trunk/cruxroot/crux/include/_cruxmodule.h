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

#include <Python.h>

#include "_cruxmodule_defs.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#ifdef HAVE_STDBOOL_H
#include <stdbool.h>
#else
#ifndef false
#define false 0
#endif
#ifndef true
#define true 1
#endif
typedef unsigned char bool;
#endif
#include <inttypes.h>
#include <unistd.h>
#include <limits.h>
#include <string.h>
#include <strings.h>
#include <sys/types.h>
#include <setjmp.h>

#include "qr.h"
#include "mem.h"
#include "xep.h"

#include "qri.h"
#include "qli.h"
#include "bhp.h"
#include "tr.h"

#ifdef CW_CPU_IA32
#include "ia32.h"
#endif

#include "tree.h"

#define CW_CRUXX_OOM 2
#define CW_CRUXX_ValueError 3

/* Used for allocation via an opaque function pointer.  These macros are used
 * to call functions such as mem_free_e(). */
#ifdef CW_DBG
#define cw_opaque_alloc(a_func, a_arg, a_size)				\
    (a_func)((void *) (a_arg), (size_t) (a_size), __FILE__, __LINE__)
#define cw_opaque_calloc(a_func, a_arg, a_num, a_size)			\
    (a_func)((void *) (a_arg), (size_t) (a_num), (size_t) (a_size),	\
	     __FILE__, __LINE__)
#define cw_opaque_realloc(a_func, a_arg, a_ptr, a_size, a_old_size)	\
    (a_func)((void *) (a_arg), (void *) (a_ptr), (size_t) (a_size),	\
	     (size_t) (a_old_size), __FILE__, __LINE__)
#define cw_opaque_dealloc(a_func, a_arg, a_ptr, a_size)			\
    (a_func)((void *) (a_arg), (void *) (a_ptr), (size_t) (a_size),	\
	     __FILE__, __LINE__)
#else
#define cw_opaque_alloc(a_func, a_arg, a_size)				\
    (a_func)((void *) (a_arg), (size_t) (a_size), NULL, 0)
#define cw_opaque_calloc(a_func, a_arg, a_num, a_size)			\
    (a_func)((void *) (a_arg), (size_t) (a_num), (size_t) (a_size),	\
	     NULL, 0)
#define cw_opaque_realloc(a_func, a_arg, a_ptr, a_size, a_old_size)	\
    (a_func)((void *) (a_arg), (void *) (a_ptr), (size_t) (a_size),	\
	     (size_t) (a_old_size), NULL, 0)
#define cw_opaque_dealloc(a_func, a_arg, a_ptr, a_size)			\
    (a_func)((void *) (a_arg), (void *) (a_ptr), (size_t) (a_size),	\
	     NULL, 0)
#endif

#ifdef WORDS_BIGENDIAN
#define cw_ntohq(a) (a)
#define cw_htonq(a) (a)
#else
#define cw_ntohq(a)							\
    (uint64_t) (((uint64_t) (ntohl((uint32_t) ((a) >> 32))))	\
		   | (((uint64_t) (ntohl((uint32_t)		\
		   ((a) & 0x00000000ffffffff)))) << 32))
#define cw_htonq(a)							\
    (uint64_t) (((uint64_t) (htonl((uint32_t) ((a) >> 32))))	\
		   | (((uint64_t) (htonl((uint32_t)		\
		   ((a) & 0x00000000ffffffff)))) << 32))
#endif

/* assert()-alike.  It's a bit prettier and cleaner, but the same idea. */
#define cw_error(a)							\
    do									\
    {									\
	fprintf(stderr, "%s:%u:%s(): Error: %s\n", __FILE__,		\
		__LINE__, __func__, a);				\
		abort();						\
    } while (0)

#ifdef CW_ASSERT
#define cw_not_reached()						\
    do									\
    {									\
	fprintf(stderr, "%s:%u:%s(): Unreachable code reached\n",	\
		__FILE__, __LINE__, __func__);			\
		abort();						\
    } while (0)

#define cw_assert(a)							\
    do									\
    {									\
	if (!(a))							\
	{								\
	    fprintf(stderr, "%s:%u:%s(): Failed assertion: \"%s\"\n",	\
		    __FILE__, __LINE__, __func__, #a);		\
	    abort();							\
	}								\
    } while (0)

/* Macro to do the drudgery of assuring that a pointer is non-NULL. */
#define cw_check_ptr(x)							\
    do									\
    {									\
	if (((x) == NULL) || ((x) == (void *) 0xa5a5a5a5)		\
	    || ((x) == (void *) 0x5a5a5a5a))				\
	{								\
	    fprintf(stderr, "%s:%u:%s(): Invalid pointer: %s (%p)\n",	\
		    __FILE__, __LINE__, __func__, #x, (x));		\
	    abort();							\
	}								\
    } while (0)
#else
#define cw_not_reached()
#define cw_assert(a)
#define cw_check_ptr(a)
#endif

/* cw_dasssert() is used internally in places that the assertion should only
 * be made if CW_DBG is defined, such as checking magic variables that only
 * exist in that case. */
#if (defined(CW_DBG) && defined(CW_ASSERT))
#define cw_dassert(a)							\
    do									\
    {									\
	if (!(a))							\
	{								\
	    fprintf(stderr, "%s:%u:%s(): Failed assertion: \"%s\"\n",	\
		    __FILE__, __LINE__, __func__, #a);		\
	    abort();							\
	}								\
    } while (0)
#else
#define cw_dassert(a)
#endif

/* Convenience macro for determining the offset of a field within a
 * structure. */
#define cw_offsetof(a_type, a_field)					\
    ((uint32_t) &(((a_type *)NULL)->a_field))
