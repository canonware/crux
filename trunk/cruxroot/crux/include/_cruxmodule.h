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

/* Python headers. */
#include <Python.h>
#include <compile.h>
#include <eval.h>

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
#include "tr.h"

#ifdef CW_CPU_IA32
#include "ia32.h"
#endif

#include "tree.h"
#include "tree_canonize.h"
#include "tree_mp.h"
#include "tree_nj.h"
#include "tree_tbr.h"

#define CW_CRUXX_OOM 2
#define CW_CRUXX_ValueError 3

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
