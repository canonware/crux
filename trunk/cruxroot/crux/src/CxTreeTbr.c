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
CxTreeTbr(CxtCxtTreeObject *self, PyObject *args)
{
    PyObject *retval
#ifdef CxmCcSilence
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

    CxmXepBegin();
    CxmXepTry
    {
	CxTrTbr(self->tr, bisect, reconnect_a, reconnect_b);

	Py_INCREF(Py_None);
	retval = Py_None;
    }
    CxmXepCatch(CxmXepOOM)
    {
	CxmXepHandled();
	retval = PyErr_NoMemory();
    }
    CxmXepEnd();

    RETURN:
    return retval;
}

PyObject *
CxTreeTbrNneighborsGet(CxtCxtTreeObject *self)
{
    PyObject *retval
#ifdef CxmCcSilence
	= NULL
#endif
	;

    CxmXepBegin();
    CxmXepTry
    {
	retval = Py_BuildValue("i", CxTrTbrNneighborsGet(self->tr));
    }
    CxmXepCatch(CxmXepOOM)
    {
	CxmXepHandled();
	retval = PyErr_NoMemory();
    }
    CxmXepEnd();

    return retval;
}

PyObject *
CxTreeTbrNeighborGet(CxtCxtTreeObject *self, PyObject *args)
{
    PyObject *retval
#ifdef CxmCcSilence
	= NULL
#endif
	;
    uint32_t neighbor, bisect, reconnect_a, reconnect_b;

    if (PyArg_ParseTuple(args, "i", &neighbor) == 0)
    {
	retval = NULL;
	goto RETURN;
    }

    CxmXepBegin();
    CxmXepTry
    {
	if (neighbor >= CxTrTbrNneighborsGet(self->tr))
	{
	    Py_INCREF(PyExc_ValueError);
	    retval = PyExc_ValueError;
	}
	else
	{
	    CxTrTbrNeighborGet(self->tr, neighbor,
				&bisect, &reconnect_a, &reconnect_b);

	    retval = Py_BuildValue("(iii)", bisect, reconnect_a, reconnect_b);
	}
    }
    CxmXepCatch(CxmXepOOM)
    {
	CxmXepHandled();
	retval = PyErr_NoMemory();
    }
    CxmXepEnd();

    RETURN:
    return retval;
}
