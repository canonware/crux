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
tree_mp_prepare(TreeObject *self, PyObject *args)
{
    PyObject *retval
#ifdef CW_CC_SILENCE
	= NULL
#endif
	;
    PyObject *taxa, *tobj;
    uint32_t elim, ntaxa, i, j;
    uint32_t nchars
#ifdef CW_CC_SILENCE
	= 0
#endif
	;
    char **tarr
#ifdef CW_CC_SILENCE
	= NULL
#endif
	;

    elim = 1;
    if (PyArg_ParseTuple(args, "O!|i", &PyList_Type, &taxa, &elim) == 0)
    {
	retval = NULL;
	goto RETURN;
    }

    ntaxa = PyList_Size(taxa);
    if (elim != 0 && elim != 1)
    {
	Py_INCREF(PyExc_ValueError);
	retval = PyExc_ValueError;
	goto RETURN;
    }

    xep_begin();
    xep_try
    {
	/* Make sure that all taxa have the same number of characters. */
	if (ntaxa > 0)
	{
	    tobj = PyList_GetItem(taxa, 0);
	    /* Don't worry about throwing a ValueError error here, since the for
	     * loop below will do so. */
	    if (PyString_Check(tobj))
	    {
		nchars = PyString_Size(tobj);
	    }

	    /* Create an array of string pointers. */
	    tarr = (char **) cw_malloc(sizeof(char *) * ntaxa);
	}

	for (i = 0; i < ntaxa; i++)
	{
	    tobj = PyList_GetItem(taxa, i);
	    if (PyString_Check(tobj) == 0 || PyString_Size(tobj) != nchars)
	    {
		cw_free(tarr);
		xep_throw(CW_CRUXX_ValueError);
	    }
	    tarr[i] = PyString_AsString(tobj);

	    /* Make sure characters are valid codes. */
	    for (j = 0; j < nchars; j++)
	    {
		switch (tarr[i][j])
		{
		    case 'N':
		    case 'n':
		    case 'X':
		    case 'x':
		    case 'V':
		    case 'v':
		    case 'H':
		    case 'h':
		    case 'M':
		    case 'm':
		    case 'D':
		    case 'd':
		    case 'R':
		    case 'r':
		    case 'W':
		    case 'w':
		    case 'A':
		    case 'a':
		    case 'B':
		    case 'b':
		    case 'S':
		    case 's':
		    case 'Y':
		    case 'y':
		    case 'C':
		    case 'c':
		    case 'K':
		    case 'k':
		    case 'G':
		    case 'g':
		    case 'T':
		    case 't':
		    case '-':
		    {
			break;
		    }
		    default:
		    {
			cw_free(tarr);
			xep_throw(CW_CRUXX_ValueError);
		    }
		}
	    }
	}

	// XXX Recurse through the tree and make sure that the taxa are numbered
	// correctly.

	/* Do preparation. */
	tr_mp_prepare(self->tr, elim, tarr, ntaxa, nchars);

	/* Clean up. */
	cw_free(tarr);

	Py_INCREF(Py_None);
	retval = Py_None;
    }
    xep_catch(CW_CRUXX_OOM)
    {
	xep_handled();
	retval = PyErr_NoMemory();
    }
    xep_catch(CW_CRUXX_ValueError)
    {
	xep_handled();
	Py_INCREF(PyExc_ValueError);
	retval = PyExc_ValueError;
    }
    xep_end();

    RETURN:
    return retval;
}

PyObject *
tree_mp_finish(TreeObject *self)
{
    tr_mp_finish(self->tr);

    Py_INCREF(Py_None);
    return Py_None;
}

PyObject *
tree_mp(TreeObject *self)
{
    return Py_BuildValue("i", tr_mp_score(self->tr));
}

PyObject *
tree_tbr_best_neighbors_mp(TreeObject *self, PyObject *args)
{
    PyObject *retval
#ifdef CW_CC_SILENCE
	= NULL
#endif
	;
    uint32_t maxhold;

    maxhold = CW_TR_HOLD_ALL;
    if (PyArg_ParseTuple(args, "|i", &maxhold) == 0)
    {
	retval = NULL;
	goto RETURN;
    }
    if (maxhold < 0 && maxhold != CW_TR_HOLD_ALL)
    {
	Py_INCREF(PyExc_ValueError);
	retval = PyExc_ValueError;
	goto RETURN;
    }

    xep_begin();
    xep_try
    {
	tr_tbr_best_neighbors_mp(self->tr, maxhold);

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
tree_tbr_better_neighbors_mp(TreeObject *self, PyObject *args)
{
    PyObject *retval
#ifdef CW_CC_SILENCE
	= NULL
#endif
	;
    uint32_t maxhold;

    maxhold = CW_TR_HOLD_ALL;
    if (PyArg_ParseTuple(args, "|i", &maxhold) == 0)
    {
	retval = NULL;
	goto RETURN;
    }
    if (maxhold < 0 && maxhold != CW_TR_HOLD_ALL)
    {
	Py_INCREF(PyExc_ValueError);
	retval = PyExc_ValueError;
	goto RETURN;
    }

    xep_begin();
    xep_try
    {
	tr_tbr_better_neighbors_mp(self->tr, maxhold);

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
tree_tbr_all_neighbors_mp(TreeObject *self)
{
    PyObject *retval
#ifdef CW_CC_SILENCE
	= NULL
#endif
	;

    xep_begin();
    xep_try
    {
	tr_tbr_all_neighbors_mp(self->tr);

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

PyObject *
tree_nheld_get(TreeObject *self)
{
    PyObject *retval
#ifdef CW_CC_SILENCE
	= NULL
#endif
	;

    xep_begin();
    xep_try
    {
	retval = Py_BuildValue("i", tr_nheld_get(self->tr));
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
tree_held_get(TreeObject *self, PyObject *args)
{
    PyObject *retval
#ifdef CW_CC_SILENCE
	= NULL
#endif
	;
    uint32_t held, neighbor, score, bisect, reconnect_a, reconnect_b;

    if (PyArg_ParseTuple(args, "i", &held) == 0)
    {
	retval = NULL;
	goto RETURN;
    }
    if (held >= tr_nheld_get(self->tr))
    {
	Py_INCREF(PyExc_ValueError);
	retval = PyExc_ValueError;
	goto RETURN;
    }

    xep_begin();
    xep_try
    {
	tr_held_get(self->tr, held, &neighbor, &score);
	tr_tbr_neighbor_get(self->tr, neighbor,
			    &bisect, &reconnect_a, &reconnect_b);

	retval = Py_BuildValue("(i(iii))", score,
			       bisect, reconnect_a, reconnect_b);
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
