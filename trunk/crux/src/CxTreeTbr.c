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
CxTreeTbr(CxtTreeObject *self, PyObject *args)
{
    PyObject *retval
#ifdef CxmCcSilence
	= NULL
#endif
	;
    uint32_t bisect, reconnectA, reconnectB;

    if (PyArg_ParseTuple(args, "(iii)", &bisect, &reconnectA, &reconnectB)
	== 0)
    {
	retval = NULL;
	goto RETURN;
    }

    CxmXepBegin();
    CxmXepTry
    {
	CxTrTbr(self->tr, bisect, reconnectA, reconnectB);

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
CxTreeTbrNneighborsGet(CxtTreeObject *self)
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
CxTreeTbrNeighborGet(CxtTreeObject *self, PyObject *args)
{
    PyObject *retval
#ifdef CxmCcSilence
	= NULL
#endif
	;
    uint32_t neighbor, bisect, reconnectA, reconnectB;

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
	    CxError(CxgTreeValueError,
		    "neighbor: %u is out of range [0..%u]",
		    neighbor, CxTrTbrNneighborsGet(self->tr));
	    retval = NULL;
	    Py_INCREF(PyExc_ValueError);
	    retval = PyExc_ValueError;
	}
	else
	{
	    CxTrTbrNeighborGet(self->tr, neighbor,
				&bisect, &reconnectA, &reconnectB);

	    retval = Py_BuildValue("(iii)", bisect, reconnectA, reconnectB);
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
