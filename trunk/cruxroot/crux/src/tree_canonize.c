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
tree_canonize(TreeObject *self)
{
    PyObject *retval
#ifdef CW_CC_SILENCE
	= NULL
#endif
	;

    xep_begin();
    xep_try
    {
	cw_tr_node_t old_tr_node, tr_node;
	NodeObject *node;

	old_tr_node = tr_base_get(self->tr);

	tr_canonize(self->tr);

	/* Reference new base. */
	tr_node = tr_base_get(self->tr);
	if (tr_node != CW_TR_NODE_NONE)
	{
	    node = (NodeObject *) tr_node_aux_get(self->tr, tr_node);
	    Py_INCREF(node);
	}

	/* Decref old base. */
	if (old_tr_node != CW_TR_NODE_NONE)
	{
	    node = (NodeObject *) tr_node_aux_get(self->tr, old_tr_node);
	    Py_DECREF(node);
	}

	Py_INCREF(Py_None);
	retval = Py_None;
    }
    xep_catch(CW_CRUXX_OOM)
    {
	xep_handled();
	retval = PyErr_NoMemory();
    }
    xep_end();

    return retval;
}
