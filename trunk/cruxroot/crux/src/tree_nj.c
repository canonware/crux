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
tree__nj(TreeObject *self, PyObject *args)
{
    PyObject *retval, *dist_list, *tobj;
    double *distances;
    uint32_t ntaxa, nelms, i;
    bool okay;

    if (PyArg_ParseTuple(args, "O!", &PyList_Type, &dist_list) == 0)
    {
	retval = NULL;
	goto RETURN;
    }

    nelms = PyList_Size(dist_list);
    ntaxa = sqrt(nelms);
    if (ntaxa * ntaxa != nelms)
    {
	Py_INCREF(PyExc_ValueError);
	retval = PyExc_ValueError;
	goto RETURN;
    }

    xep_begin();
    xep_try
    {
	distances = (double *) cw_malloc(sizeof(double) * nelms);

	okay = true;
	for (i = 0; i < nelms; i++)
	{
	    tobj = PyList_GetItem(dist_list, i);
	    if (PyFloat_Check(tobj))
	    {
		distances[i] = PyFloat_AsDouble(tobj);
	    }
	    else if (PyInt_Check(tobj))
	    {
		distances[i] = PyInt_AsLong(tobj);
	    }
	    else
	    {
		Py_INCREF(PyExc_ValueError);
		retval = PyExc_ValueError;
		okay = false;
	    }
	}

	if (okay)
	{
	    cw_tr_node_t old_tr_node, tr_node;
	    NodeObject *node;

	    old_tr_node = tr_base_get(self->tr);

	    /* Neighbor-join. */
	    tr_nj(self->tr, distances, ntaxa);
	    
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
	}

	cw_free(distances);
    }
    xep_catch(CW_CRUXX_OOM)
    {
	xep_handled();
	retval = PyErr_NoMemory();
    }
    xep_end();

    Py_INCREF(Py_None);
    retval = Py_None;
    RETURN:
    return retval;
}
