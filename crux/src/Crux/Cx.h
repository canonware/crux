#ifndef Cx_h
#define Cx_h

// XXX Remove Python header includes?
// Python headers.
//#include <Python.h>
//#include <compile.h>
//#include <eval.h>

#include "Cx_defs.h"

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

// assert()-alike.  It's a bit prettier and cleaner, but the same idea.
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

// Macro to do the drudgery of assuring that a pointer is non-NULL.
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

// CxmDassert() is used internally in places that the assertion should only be
// made if CxmDebug is defined, such as checking magic variables that only exist
// in that case.
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

// Convenience macro for determining the offset of a field within a structure.
#define CxmOffsetOf(aType, aField)					\
    ((size_t) &(((aType *)NULL)->aField))

#include "CxIa32.h"
#include "CxAmd64.h"
#include "CxPpc.h"

void
CxInit(void);

#endif // Cx_h
