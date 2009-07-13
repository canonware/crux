#ifndef Cx_h
#define Cx_h

#include <Cx_defs.h>

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
#include <pthread.h>

#ifdef CxmMpi
#  include <mpi.h>
#endif

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
	if ((x) == NULL)						\
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

extern unsigned CxNcpus;

void
CxInit(void);
void
CxThreaded(void);

#ifndef CxmUseInlines
int
CxCmp2Richcmp(int cmp, int op);
#endif

#if (defined(CxmUseInlines) || defined(Cx_c))
// Convert {-1,0,1} for the __richcmp__ special method.
CxmInline int
CxCmp2Richcmp(int cmp, int op) {
    static const int tab[] = {
	/*
	-1, 0, 1, x */
	1, 0, 0, -1, // <
	1, 1, 0, -1, // <=
	0, 1, 0, -1, // ==
	1, 0, 1, -1, // !=
	0, 0, 1, -1, // >
	0, 1, 1, -1  // >=
    };

    CxmAssert(cmp >= -1 && cmp <= 1);
    CxmAssert(op >= 0 && op < 6);

    return tab[(op << 2) + (cmp + 1)];
}
#endif

#endif // Cx_h
