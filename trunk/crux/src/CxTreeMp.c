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
CxTreeMpPrepare(CxtTreeObject *self, PyObject *args)
{
    PyObject *retval
#ifdef CxmCcSilence
	= NULL
#endif
	;
    PyObject *taxa, *tobj;
    uint32_t elim, ntaxa, i, j;
    uint32_t nchars
#ifdef CxmCcSilence
	= 0
#endif
	;
    char **tarr
#ifdef CxmCcSilence
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
	CxError(CxgTreeValueError,
		"elim: False or True expected");
	retval = NULL;
	goto RETURN;
    }

    CxmXepBegin();
    CxmXepTry
    {
	/* Make sure that all taxa have the same number of characters. */
	if (ntaxa > 0)
	{
	    tobj = PyList_GetItem(taxa, 0);
	    /* Don't worry about raising ValueError error here, since the for
	     * loop below will do so. */
	    if (PyString_Check(tobj))
	    {
		nchars = PyString_Size(tobj);
	    }

	    /* Create an array of string pointers. */
	    tarr = (char **) CxmMalloc(sizeof(char *) * ntaxa);
	}

	for (i = 0; i < ntaxa; i++)
	{
	    tobj = PyList_GetItem(taxa, i);
	    if (PyString_Check(tobj) == 0 || PyString_Size(tobj) != nchars)
	    {
		CxError(CxgTreeValueError,
			"Character string expected");
		goto ERROR;
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
			CxError(CxgTreeValueError,
				"Invalid character '%c'", tarr[i][j]);
			goto ERROR;
		    }
		}
	    }
	}

	// XXX Recurse through the tree and make sure that the taxa are numbered
	// correctly.

	/* Do preparation. */
	CxTrMpPrepare(self->tr, elim, tarr, ntaxa, nchars);

	/* Clean up. */
	CxmFree(tarr);

	Py_INCREF(Py_None);
	retval = Py_None;
	break;

	ERROR:
	CxmFree(tarr);
	retval = NULL;
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
CxTreeMpFinish(CxtTreeObject *self)
{
    CxTrMpFinish(self->tr);

    Py_INCREF(Py_None);
    return Py_None;
}

PyObject *
CxTreeMp(CxtTreeObject *self)
{
    return Py_BuildValue("i", CxTrMpScore(self->tr));
}

PyObject *
CxTreeTbrBestNeighborsMp(CxtTreeObject *self, PyObject *args)
{
    PyObject *retval
#ifdef CxmCcSilence
	= NULL
#endif
	;
    int maxHold;

    maxHold = CxmTrHoldAll;
    if (PyArg_ParseTuple(args, "|i", &maxHold) == 0)
    {
	retval = NULL;
	goto RETURN;
    }
    if (maxHold < 0 && maxHold != CxmTrHoldAll)
    {
	CxError(CxgTreeValueError,
		"maxHold: non-negative integer expected");
	retval = NULL;
	goto RETURN;
    }

    CxmXepBegin();
    CxmXepTry
    {
	CxTrTbrBestNeighborsMp(self->tr, maxHold);

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
CxTreeTbrBetterNeighborsMp(CxtTreeObject *self, PyObject *args)
{
    PyObject *retval
#ifdef CxmCcSilence
	= NULL
#endif
	;
    int maxHold;

    maxHold = CxmTrHoldAll;
    if (PyArg_ParseTuple(args, "|i", &maxHold) == 0)
    {
	retval = NULL;
	goto RETURN;
    }
    if (maxHold < 0 && maxHold != CxmTrHoldAll)
    {
	CxError(CxgTreeValueError,
		"maxHold: non-negative integer expected");
	retval = NULL;
	goto RETURN;
    }

    CxmXepBegin();
    CxmXepTry
    {
	CxTrTbrBetterNeighborsMp(self->tr, maxHold);

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
CxTreeTbrAllNeighborsMp(CxtTreeObject *self)
{
    PyObject *retval
#ifdef CxmCcSilence
	= NULL
#endif
	;

    CxmXepBegin();
    CxmXepTry
    {
	CxTrTbrAllNeighborsMp(self->tr);

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

PyObject *
CxTreeNheldGet(CxtTreeObject *self)
{
    PyObject *retval
#ifdef CxmCcSilence
	= NULL
#endif
	;

    CxmXepBegin();
    CxmXepTry
    {
	retval = Py_BuildValue("i", CxTrNheldGet(self->tr));
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
CxTreeheldGet(CxtTreeObject *self, PyObject *args)
{
    PyObject *retval
#ifdef CxmCcSilence
	= NULL
#endif
	;
    int held;
    uint32_t neighbor, score, bisect, reconnectA, reconnectB;

    if (PyArg_ParseTuple(args, "i", &held) == 0)
    {
	retval = NULL;
	goto RETURN;
    }
    if (held >= CxTrNheldGet(self->tr))
    {
	Py_INCREF(PyExc_ValueError);
	retval = PyExc_ValueError;
	goto RETURN;
    }

    CxmXepBegin();
    CxmXepTry
    {
	CxTrHeldGet(self->tr, held, &neighbor, &score);
	CxTrTbrNeighborGet(self->tr, neighbor,
			    &bisect, &reconnectA, &reconnectB);

	retval = Py_BuildValue("(i(iii))", score,
			       bisect, reconnectA, reconnectB);
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
