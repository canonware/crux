//==============================================================================
//
// <Copyright = jasone>
// <License>
//
//==============================================================================
//
// Version: Crux <Version = crux>
//
//==============================================================================

#include "../include/_cruxmodule.h"

PyObject *
CxTreeTbr(CxtTreeObject *self, PyObject *args)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    uint32_t bisect, reconnectA, reconnectB;

    if (PyArg_ParseTuple(args, "(iii)", &bisect, &reconnectA, &reconnectB)
	== 0)
    {
	rVal = NULL;
	goto RETURN;
    }

    CxmXepBegin();
    CxmXepTry
    {
	CxTrTbr(self->tr, bisect, reconnectA, reconnectB);

	Py_INCREF(Py_None);
	rVal = Py_None;
    }
    CxmXepCatch(CxmXepOOM)
    {
	CxmXepHandled();
	rVal = PyErr_NoMemory();
    }
    CxmXepEnd();

    RETURN:
    return rVal;
}

PyObject *
CxTreeTbrNneighborsGet(CxtTreeObject *self)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;

    CxmXepBegin();
    CxmXepTry
    {
	rVal = Py_BuildValue("i", CxTrTbrNneighborsGet(self->tr));
    }
    CxmXepCatch(CxmXepOOM)
    {
	CxmXepHandled();
	rVal = PyErr_NoMemory();
    }
    CxmXepEnd();

    return rVal;
}

PyObject *
CxTreeTbrNeighborGet(CxtTreeObject *self, PyObject *args)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    uint32_t neighbor, bisect, reconnectA, reconnectB;

    if (PyArg_ParseTuple(args, "i", &neighbor) == 0)
    {
	rVal = NULL;
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
	    rVal = NULL;
	    Py_INCREF(PyExc_ValueError);
	    rVal = PyExc_ValueError;
	}
	else
	{
	    CxTrTbrNeighborGet(self->tr, neighbor,
				&bisect, &reconnectA, &reconnectB);

	    rVal = Py_BuildValue("(iii)", bisect, reconnectA, reconnectB);
	}
    }
    CxmXepCatch(CxmXepOOM)
    {
	CxmXepHandled();
	rVal = PyErr_NoMemory();
    }
    CxmXepEnd();

    RETURN:
    return rVal;
}
