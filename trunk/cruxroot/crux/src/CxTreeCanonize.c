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
CxTreeCanonize(CxtCxtTreeObject *self)
{
    PyObject *retval
#ifdef CxmCcSilence
	= NULL
#endif
	;

    CxmXepBegin();
    CxmXepTry
    {
	CxtTrNode old_tr_node, tr_node;
	CxtCxtNodeObject *node;

	old_tr_node = CxTrBaseGet(self->tr);

	CxTrCanonize(self->tr);

	/* Reference new base. */
	tr_node = CxTrBaseGet(self->tr);
	if (tr_node != CxmTrNodeNone)
	{
	    node = (CxtCxtNodeObject *) CxTrNodeAuxGet(self->tr, tr_node);
	    Py_INCREF(node);
	}

	/* Decref old base. */
	if (old_tr_node != CxmTrNodeNone)
	{
	    node = (CxtCxtNodeObject *) CxTrNodeAuxGet(self->tr, old_tr_node);
	    Py_DECREF(node);
	}

	Py_INCREF(Py_None);
	retval = Py_None;
    }
    CxmXepCatch(CxmXepOOM)
    {
	CxmXepHandled();
	retval = PyErr_NoMemory();
    }
    CxmXepEnd();

    return retval;
}
