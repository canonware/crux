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

#ifdef CxmDebug
static bool CxpXepInitialized = false;
#endif
static CxtXep *CxpXepFirst;

void
CxXepInit(void)
{
    cxmAssert(CxpXepInitialized == false);

    CxpXepFirst = NULL;
#ifdef CxmDebug
    CxpXepInitialized = true;
#endif
}

void
CxXepShutdown(void)
{
    cxmAssert(CxpXepInitialized);

#ifdef CxmDebug
    CxpXepInitialized = false;
#endif
}

void
CxXepThrowE(CxtXepv aValue, volatile const char *aFilename,
	    uint32_t aLineNum)
{
    CxtXep *xepFirst, *xep;

    cxmAssert(CxpXepInitialized);
    cxmAssert(aValue > CxeXepsCatch);

    /* Iterate backward through the exception handlers until the exception is
     * handled or there are no more exception handlers. */
    xep = xepFirst = CxpXepFirst;
    if (xepFirst != NULL)
    {
	xep = CxmQrPrev(xepFirst, link);
    }
    else
    {
	/* No exception handlers at all. */
	fprintf(stderr, "%s(): Unhandled exception %u thrown at %s:%u\n",
		__func__, aValue, aFilename, aLineNum);
	abort();
    }

    do
    {
	xep->isHandled = false;
	xep->filename = aFilename;
	xep->lineNum = aLineNum;

	switch (xep->state)
	{
	    case CxeXepsTry:
	    {
		/* Execute the handler. */
		xep->value = aValue;
		xep->state = CxeXepsCatch;
		longjmp(xep->context, (int) aValue);
		cxmNotReached();
	    }
	    case CxeXepsCatch:
	    {
		/* Exception thrown within handler; propagate. */
		break;
	    }
	    default:
	    {
		cxmNotReached();
	    }
	}

	xep = CxmQrPrev(xep, link);
    } while (xep != xepFirst);

    /* No more exception handlers. */
    fprintf(stderr, "%s(): Unhandled exception %u thrown at %s:%u\n",
	    __func__, aValue, xep->filename, xep->lineNum);
    abort();
}

void
CxpXepRetry(CxtXep *aXep)
{
    cxmAssert(CxpXepInitialized);

#ifdef CxmDebug
    switch (aXep->state)
    {
	case CxeXepsCatch:
	{
	    break;
	}
	case CxeXepsTry:
	{
	    cxmError("Exception retry outside handler");
	}
	default:
	{
	    cxmNotReached();
	}
    }
#endif
    aXep->value = CxmXepvNone;
    aXep->state = CxeXepsTry;
    aXep->isHandled = true;
    longjmp(aXep->context, (int) CxmXepvCode);
    cxmNotReached();
}

void
CxpXepHandled(CxtXep *aXep)
{
    cxmAssert(CxpXepInitialized);

#ifdef CxmDebug
    switch (aXep->state)
    {
	case CxeXepsCatch:
	{
	    break;
	}
	case CxeXepsTry:
	{
	    cxmError("Exception handled outside handler");
	}
	default:
	{
	    cxmNotReached();
	}
    }
#endif

    aXep->isHandled = true;
    CxpXepUnlink(aXep);
}

void
CxpXepLink(CxtXep *aXep)
{
    CxtXep *xepFirst;

    cxmAssert(CxpXepInitialized);

    xepFirst = CxpXepFirst;

    /* Link into the xep ring, if it exists. */
    CxmQrNew(aXep, link);
    if (xepFirst != NULL)
    {
	cxmCheckPtr(CxmQrPrev(xepFirst, link));
	cxmCheckPtr(CxmQrNext(xepFirst, link));

	CxmQrBeforeInsert(xepFirst, aXep, link);
    }
    else
    {
	CxpXepFirst = aXep;
    }

    aXep->value = CxmXepvNone;
    aXep->state = CxeXepsTry;
    aXep->isHandled = true;
    aXep->isLinked = true;
}

void
CxpXepUnlink(CxtXep *aXep)
{
    CxtXep *xepFirst;

    cxmAssert(CxpXepInitialized);

    if (aXep->isLinked)
    {
	xepFirst = CxpXepFirst;
	cxmCheckPtr(CxmQrPrev(xepFirst, link));
	cxmCheckPtr(CxmQrNext(xepFirst, link));

	/* Remove handler from ring. */
	if (aXep != xepFirst)
	{
	    CxmQrRemove(aXep, link);
	}
	else
	{
	    CxpXepFirst = NULL;
	}
	aXep->isLinked = false;

	if (aXep->isHandled == false)
	{
	    if (aXep != xepFirst)
	    {
		/* Propagate exception. */
		CxXepThrowE(aXep->value, aXep->filename,
			    aXep->lineNum);
	    }
	    else
	    {
		/* No more exception handlers. */
		fprintf(stderr, "%s(): Unhandled exception "
			"%u thrown at %s:%u\n", __func__,
			aXep->value, aXep->filename,
			aXep->lineNum);
		abort();
	    }
	}
    }
}
