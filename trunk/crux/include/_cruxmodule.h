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

#define CxmXepOOM 2
#define CxmXepValueError 3

extern PyObject *CxgException;

/* assert()-alike.  It's a bit prettier and cleaner, but the same idea. */
#define CxmError(a)							\
    do									\
    {									\
	fprintf(stderr, "%s:%d:%s(): Error: %s\n", __FILE__,		\
		__LINE__, __func__, a);					\
		abort();						\
    } while (0)

#ifdef CxmAssertions
#define CxmNotReached()							\
    do									\
    {									\
	fprintf(stderr, "%s:%d:%s(): Unreachable code reached\n",	\
		__FILE__, __LINE__, __func__);				\
		abort();						\
    } while (0)

#define CxmAssert(a)							\
    do									\
    {									\
	if (!(a))							\
	{								\
	    fprintf(stderr, "%s:%d:%s(): Failed assertion: \"%s\"\n",	\
		    __FILE__, __LINE__, __func__, #a);			\
	    abort();							\
	}								\
    } while (0)

/* Macro to do the drudgery of assuring that a pointer is non-NULL. */
#define CxmCheckPtr(x)							\
    do									\
    {									\
	if (((x) == NULL) || ((x) == (void *) 0xa5a5a5a5)		\
	    || ((x) == (void *) 0x5a5a5a5a))				\
	{								\
	    fprintf(stderr, "%s:%d:%s(): Invalid pointer: %s (%p)\n",	\
		    __FILE__, __LINE__, __func__, #x, (x));		\
	    abort();							\
	}								\
    } while (0)
#else
#define CxmNotReached()
#define CxmAssert(a)
#define CxmCheckPtr(a)
#endif

/* CxmDassert() is used internally in places that the assertion should only
 * be made if CxmDebug is defined, such as checking magic variables that only
 * exist in that case. */
#if (defined(CxmDebug) && defined(CxmAssertions))
#define CxmDassert(a)							\
    do									\
    {									\
	if (!(a))							\
	{								\
	    fprintf(stderr, "%s:%d:%s(): Failed assertion: \"%s\"\n",	\
		    __FILE__, __LINE__, __func__, #a);			\
	    abort();							\
	}								\
    } while (0)
#else
#define CxmDassert(a)
#endif

/* Convenience macro for determining the offset of a field within a
 * structure. */
#define CxmOffsetOf(aType, aField)					\
    ((uint32_t) &(((aType *)NULL)->aField))

#include "CxMt.h"
#include "CxRi.h"
#include "CxQr.h"
#include "CxMem.h"
#include "CxXep.h"

#include "CxQri.h"
#include "CxQli.h"
#include "CxTr.h"

#ifdef CxmCpuIa32
#include "CxIa32.h"
#endif

#include "CxTree.h"
#include "CxTreeCanonize.h"
#include "CxTreeMp.h"
#include "CxTreeNj.h"
#include "CxTreeTbr.h"
#include "CxFastaParser.h"
#include "CxDistMatrix.h"

void
init_crux(void);

void
CxError(PyObject *exception, const char *format, ...);
