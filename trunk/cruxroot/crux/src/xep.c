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

#ifdef CW_DBG
static bool s_xep_initialized = false;
#endif
static cw_xep_t *s_xep_first;

void
xep_init(void)
{
    cw_assert(s_xep_initialized == false);

    s_xep_first = NULL;
#ifdef CW_DBG
    s_xep_initialized = true;
#endif
}

void
xep_shutdown(void)
{
    cw_assert(s_xep_initialized);

#ifdef CW_DBG
    s_xep_initialized = false;
#endif
}

void
xep_throw_e(cw_xepv_t a_value, volatile const char *a_filename,
	    uint32_t a_line_num)
{
    cw_xep_t *xep_first, *xep;

    cw_assert(s_xep_initialized);
    cw_assert(a_value > CW_XEPS_CATCH);

    /* Iterate backward through the exception handlers until the exception is
     * handled or there are no more exception handlers. */
    xep = xep_first = s_xep_first;
    if (xep_first != NULL)
    {
	xep = qr_prev(xep_first, link);
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
	    case CW_XEPS_TRY:
	    {
		/* Execute the handler. */
		xep->value = a_value;
		xep->state = CW_XEPS_CATCH;
		longjmp(xep->context, (int) a_value);
		cw_not_reached();
	    }
	    case CW_XEPS_CATCH:
	    {
		/* Exception thrown within handler; propagate. */
		break;
	    }
	    default:
	    {
		cw_not_reached();
	    }
	}

	xep = qr_prev(xep, link);
    } while (xep != xep_first);

    /* No more exception handlers. */
    fprintf(stderr, "%s(): Unhandled exception %u thrown at %s:%u\n",
	    __func__, a_value, xep->filename, xep->line_num);
    abort();
}

void
xep_p_retry(cw_xep_t *a_xep)
{
    cw_assert(s_xep_initialized);

#ifdef CW_DBG
    switch (a_xep->state)
    {
	case CW_XEPS_CATCH:
	{
	    break;
	}
	case CW_XEPS_TRY:
	{
	    cw_error("Exception retry outside handler");
	}
	default:
	{
	    cw_not_reached();
	}
    }
#endif
    a_xep->value = CW_XEPV_NONE;
    a_xep->state = CW_XEPS_TRY;
    a_xep->is_handled = true;
    longjmp(a_xep->context, (int) CW_XEPV_CODE);
    cw_not_reached();
}

void
xep_p_handled(cw_xep_t *a_xep)
{
    cw_assert(s_xep_initialized);

#ifdef CW_DBG
    switch (a_xep->state)
    {
	case CW_XEPS_CATCH:
	{
	    break;
	}
	case CW_XEPS_TRY:
	{
	    cw_error("Exception handled outside handler");
	}
	default:
	{
	    cw_not_reached();
	}
    }
#endif

    a_xep->is_handled = true;
    xep_p_unlink(a_xep);
}

void
xep_p_link(cw_xep_t *a_xep)
{
    cw_xep_t *xep_first;

    cw_assert(s_xep_initialized);

    xep_first = s_xep_first;

    /* Link into the xep ring, if it exists. */
    qr_new(a_xep, link);
    if (xep_first != NULL)
    {
	cw_check_ptr(qr_prev(xep_first, link));
	cw_check_ptr(qr_next(xep_first, link));

	qr_before_insert(xep_first, a_xep, link);
    }
    else
    {
	s_xep_first = a_xep;
    }

    a_xep->value = CW_XEPV_NONE;
    a_xep->state = CW_XEPS_TRY;
    a_xep->is_handled = true;
    a_xep->is_linked = true;
}

void
xep_p_unlink(cw_xep_t *a_xep)
{
    cw_xep_t *xep_first;

    cw_assert(s_xep_initialized);

    if (a_xep->is_linked)
    {
	xep_first = s_xep_first;
	cw_check_ptr(qr_prev(xep_first, link));
	cw_check_ptr(qr_next(xep_first, link));

	/* Remove handler from ring. */
	if (a_xep != xep_first)
	{
	    qr_remove(a_xep, link);
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
		xep_throw_e(a_xep->value, a_xep->filename,
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
