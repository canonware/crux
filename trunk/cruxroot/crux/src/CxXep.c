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
static bool s_CxXepInitialized = false;
#endif
static CxtXep *s_xep_first;

void
CxXepInit(void)
{
    cxmAssert(s_CxXepInitialized == false);

    s_xep_first = NULL;
#ifdef CxmDebug
    s_CxXepInitialized = true;
#endif
}

void
CxXepShutdown(void)
{
    cxmAssert(s_CxXepInitialized);

#ifdef CxmDebug
    s_CxXepInitialized = false;
#endif
}

void
CxXepThrowE(CxtXepv a_value, volatile const char *a_filename,
	    uint32_t a_line_num)
{
    CxtXep *xep_first, *xep;

    cxmAssert(s_CxXepInitialized);
    cxmAssert(a_value > CxeXepsCatch);

    /* Iterate backward through the exception handlers until the exception is
     * handled or there are no more exception handlers. */
    xep = xep_first = s_xep_first;
    if (xep_first != NULL)
    {
	xep = CxmQrPrev(xep_first, link);
    }
    else
    {
	/* No exception handlers at all. */
	fprintf(stderr, "%s(): Unhandled exception %u thrown at %s:%u\n",
		__func__, a_value, a_filename, a_line_num);
	abort();
    }

    do
    {
	xep->is_handled = false;
	xep->filename = a_filename;
	xep->line_num = a_line_num;

	switch (xep->state)
	{
	    case CxeXepsTry:
	    {
		/* Execute the handler. */
		xep->value = a_value;
		xep->state = CxeXepsCatch;
		longjmp(xep->context, (int) a_value);
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
    } while (xep != xep_first);

    /* No more exception handlers. */
    fprintf(stderr, "%s(): Unhandled exception %u thrown at %s:%u\n",
	    __func__, a_value, xep->filename, xep->line_num);
    abort();
}

void
CxpXepRetry(CxtXep *a_xep)
{
    cxmAssert(s_CxXepInitialized);

#ifdef CxmDebug
    switch (a_xep->state)
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
    a_xep->value = CxmXepvNone;
    a_xep->state = CxeXepsTry;
    a_xep->is_handled = true;
    longjmp(a_xep->context, (int) CxmXepvCode);
    cxmNotReached();
}

void
CxpXepHandled(CxtXep *a_xep)
{
    cxmAssert(s_CxXepInitialized);

#ifdef CxmDebug
    switch (a_xep->state)
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

    a_xep->is_handled = true;
    CxpXepUnlink(a_xep);
}

void
CxpXepLink(CxtXep *a_xep)
{
    CxtXep *xep_first;

    cxmAssert(s_CxXepInitialized);

    xep_first = s_xep_first;

    /* Link into the xep ring, if it exists. */
    CxmQrNew(a_xep, link);
    if (xep_first != NULL)
    {
	cxmCheckPtr(CxmQrPrev(xep_first, link));
	cxmCheckPtr(CxmQrNext(xep_first, link));

	CxmQrBeforeInsert(xep_first, a_xep, link);
    }
    else
    {
	s_xep_first = a_xep;
    }

    a_xep->value = CxmXepvNone;
    a_xep->state = CxeXepsTry;
    a_xep->is_handled = true;
    a_xep->is_linked = true;
}

void
CxpXepUnlink(CxtXep *a_xep)
{
    CxtXep *xep_first;

    cxmAssert(s_CxXepInitialized);

    if (a_xep->is_linked)
    {
	xep_first = s_xep_first;
	cxmCheckPtr(CxmQrPrev(xep_first, link));
	cxmCheckPtr(CxmQrNext(xep_first, link));

	/* Remove handler from ring. */
	if (a_xep != xep_first)
	{
	    CxmQrRemove(a_xep, link);
	}
	else
	{
	    s_xep_first = NULL;
	}
	a_xep->is_linked = false;

	if (a_xep->is_handled == false)
	{
	    if (a_xep != xep_first)
	    {
		/* Propagate exception. */
		CxXepThrowE(a_xep->value, a_xep->filename,
			    a_xep->line_num);
	    }
	    else
	    {
		/* No more exception handlers. */
		fprintf(stderr, "%s(): Unhandled exception "
			"%u thrown at %s:%u\n", __func__,
			a_xep->value, a_xep->filename,
			a_xep->line_num);
		abort();
	    }
	}
    }
}
