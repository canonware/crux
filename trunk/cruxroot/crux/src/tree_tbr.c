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

PyObject *
tree_tbr(TreeObject *self, PyObject *args)
{
    PyObject *retval
#ifdef CW_CC_SILENCE
	= NULL
#endif
	;
    uint32_t bisect, reconnect_a, reconnect_b;

    if (PyArg_ParseTuple(args, "(iii)", &bisect, &reconnect_a, &reconnect_b)
	== 0)
    {
	retval = NULL;
	goto RETURN;
    }

    xep_begin();
    xep_try
    {
	tr_tbr(self->tr, bisect, reconnect_a, reconnect_b);

	Py_INCREF(Py_None);
	retval = Py_None;
    }
    xep_catch(CW_CRUXX_OOM)
    {
	xep_handled();
	retval = PyErr_NoMemory();
    }
    xep_end();

    RETURN:
    return retval;
}

PyObject *
tree_tbr_nneighbors_get(TreeObject *self)
{
    PyObject *retval
#ifdef CW_CC_SILENCE
	= NULL
#endif
	;

    xep_begin();
    xep_try
    {
	retval = Py_BuildValue("i", tr_tbr_nneighbors_get(self->tr));
    }
    xep_catch(CW_CRUXX_OOM)
    {
	xep_handled();
	retval = PyErr_NoMemory();
    }
    xep_end();

    return retval;
}

PyObject *
tree_tbr_neighbor_get(TreeObject *self, PyObject *args)
{
    PyObject *retval
#ifdef CW_CC_SILENCE
	= NULL
#endif
	;
    uint32_t neighbor, bisect, reconnect_a, reconnect_b;

    if (PyArg_ParseTuple(args, "i", &neighbor) == 0)
    {
	retval = NULL;
	goto RETURN;
    }

    xep_begin();
    xep_try
    {
	if (neighbor >= tr_tbr_nneighbors_get(self->tr))
	{
	    Py_INCREF(PyExc_ValueError);
	    retval = PyExc_ValueError;
	}
	else
	{
	    tr_tbr_neighbor_get(self->tr, neighbor,
				&bisect, &reconnect_a, &reconnect_b);

	    retval = Py_BuildValue("(iii)", bisect, reconnect_a, reconnect_b);
	}
    }
    xep_catch(CW_CRUXX_OOM)
    {
	xep_handled();
	retval = PyErr_NoMemory();
    }
    xep_end();

    RETURN:
    return retval;
}
