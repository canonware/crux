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

typedef struct CxsTreeMpData CxtTreeMpData;
typedef struct CxsTreeMpPs CxtTreeMpPs;

struct CxsTreeMpData
{
    PyObject *cTMatrix;
    uint64_t cTMatrixSeq;

    bool eliminateUninformative;

    unsigned edgeAuxInd;
    unsigned ringAuxInd;

    // XXX
};

// Character (in the systematics sense of the word).
typedef char CxtTreeMpC;

// Partial parsimony score information.
struct CxsTreeMpPs
{
    // Parent which most recently used this node's partial score when caching
    // its results.  Both children must still point to the parent in order for
    // the cached results to be valid.
    CxtTreeMpPs *parent;

    // Sum of the subtree scores, and this node's score, given particular
    // children.  In order for this to be useful, both childrens' parent
    // pointers must still point to this node.
    unsigned subtreesScore;
    unsigned nodeScore;

    // chars points to an array of Fitch parsimony state sets.  Each element in
    // the array contains a bitmap representation of a subset of {ACGT} in the 4
    // least significant bits.  T is the least significant bit.  1 means that a
    // nucleotide is in the set.
    //
    // There are nchars character state sets.
    //
    // achars is the actual allocation, which is padded in order to
    // be able to guarantee that chars is 16 byte-aligned.
    CxtTreeMpC *chars;
    unsigned nChars;
    CxtTreeMpC *aChars;
};

static void
CxpTreeMpCleanupFinal(CxtTreeObject *aTree, void *aData, unsigned aInd)
{
    CxtTreeMpData *data = (CxtTreeMpData *) aData;

    // XXX

    free(data);
}

static bool
CxpTreeMpInitEdge(CxtEdgeObject *aEdge, unsigned aInd)
{
    bool rVal;

    CxmError("XXX Not implemented");

    return rVal;
}

static void
CxpTreeMpCleanupEdge(CxtEdgeObject *aEdge, void *aData, unsigned aInd)
{
    CxmError("XXX Not implemented");
}

static bool
CxpTreeMpInitRing(CxtRingObject *aRing, unsigned aInd)
{
    bool rVal;

    CxmError("XXX Not implemented");

    return rVal;
}

static void
CxpTreeMpCleanupRing(CxtRingObject *aRing, void *aData, unsigned aInd)
{
    CxmError("XXX Not implemented");
}

static bool
CxpTreeMpDataGet(CxtTreeObject *self, CxtTreeMpData **rData)
{
    bool rVal;
    unsigned treeAuxInd;
    CxtTreeMpData *data;

    // Get aux indices for MP data.
    if (CxTreeAuxSearch(self, "MP", &treeAuxInd))
    {
	// No aux registration.
	data = (CxtTreeMpData *) malloc(sizeof(CxtTreeMpData));
	if (data == NULL)
	{
	    rVal = true;
	    goto RETURN;
	}

	// Initialize data.
	data->cTMatrix = NULL;
	data->cTMatrixSeq = 0;
	data->eliminateUninformative = false;

	// Create MP-specific aux mappings for edges and rings.
	if (CxTreeEdgeAuxRegister(self, "MP", NULL,
				  CxpTreeMpInitEdge, NULL,
				  CxpTreeMpCleanupEdge,
				  &data->edgeAuxInd)
	    || CxTreeRingAuxRegister(self, "MP", NULL,
				     CxpTreeMpInitRing, NULL,
				     CxpTreeMpCleanupRing,
				     &data->ringAuxInd)
	    || CxTreeAuxRegister(self, "MP", (void *) data,
				 CxpTreeMpCleanupFinal, NULL, &treeAuxInd))
	{
	    free(data);
	    rVal = true;
	    goto RETURN;
	}
    }
    else
    {
	data = (CxtTreeMpData *) CxTreeAuxData(self, treeAuxInd);
    }

    *rData = data;
    rVal = false;
    RETURN:
    return rVal;
}



PyObject *
CxTreeMpPrepare(CxtTreeObject *self, PyObject *args)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    PyObject *taxa, *tobj;
    uint32_t elim, ntaxa, i, j;
    char **tarr
#ifdef CxmCcSilence
	= NULL
#endif
	;

    if (PyArg_ParseTuple(args, "O!i", &PyList_Type, &taxa, &elim) == 0)
    {
	rVal = NULL;
	goto RETURN;
    }

    ntaxa = PyList_Size(taxa);

    // Make sure that all taxa have the same number of characters.
    if (ntaxa > 0)
    {
	uint32_t nchars;

	tobj = PyList_GetItem(taxa, 0);
	// Don't worry about raising ValueError error here, since the for
	// loop below will do so.
	if (PyString_Check(tobj))
	{
	    nchars = PyString_Size(tobj);
	}

	// Create an array of string pointers.
	tarr = (char **) malloc(sizeof(char *) * ntaxa);
	if (tarr == NULL)
	{
	    rVal = PyErr_NoMemory();
	    goto RETURN;
	}

	for (i = 0; i < ntaxa; i++)
	{
	    tobj = PyList_GetItem(taxa, i);
	    if (PyString_Check(tobj) == 0 || PyString_Size(tobj) != nchars)
	    {
		free(tarr);
		CxError(CxgTreeValueError,
			"Character string expected");
		rVal = NULL;
		goto RETURN;
	    }
	    tarr[i] = PyString_AsString(tobj);

	    // Make sure characters are valid codes.
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
			free(tarr);
			CxError(CxgTreeValueError,
				"Invalid character '%c'", tarr[i][j]);
			rVal = NULL;
			goto RETURN;
		    }
		}
	    }
	}

	// XXX Recurse through the tree and make sure that the taxa are
	// numbered correctly.

	// Do preparation.
	CxTrMpPrepare(self->tr, elim, tarr, ntaxa, nchars);

	// Clean up.
	free(tarr);
    }

    Py_INCREF(Py_None);
    rVal = Py_None;
    RETURN:
    return rVal;
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
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    int maxHold;

    maxHold = CxmTrHoldAll;
    if (PyArg_ParseTuple(args, "|i", &maxHold) == 0)
    {
	rVal = NULL;
	goto RETURN;
    }
    if (maxHold < 0 && maxHold != CxmTrHoldAll)
    {
	CxError(CxgTreeValueError,
		"maxHold: non-negative integer expected");
	rVal = NULL;
	goto RETURN;
    }

    CxmXepBegin();
    CxmXepTry
    {
	CxTrTbrBestNeighborsMp(self->tr, maxHold);

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
CxTreeTbrBetterNeighborsMp(CxtTreeObject *self, PyObject *args)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    int maxHold;

    maxHold = CxmTrHoldAll;
    if (PyArg_ParseTuple(args, "|i", &maxHold) == 0)
    {
	rVal = NULL;
	goto RETURN;
    }
    if (maxHold < 0 && maxHold != CxmTrHoldAll)
    {
	CxError(CxgTreeValueError,
		"maxHold: non-negative integer expected");
	rVal = NULL;
	goto RETURN;
    }

    CxmXepBegin();
    CxmXepTry
    {
	CxTrTbrBetterNeighborsMp(self->tr, maxHold);

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
CxTreeTbrAllNeighborsMp(CxtTreeObject *self)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;

    CxmXepBegin();
    CxmXepTry
    {
	CxTrTbrAllNeighborsMp(self->tr);

	Py_INCREF(Py_None);
	rVal = Py_None;
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
CxTreeNheldGet(CxtTreeObject *self)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;

    CxmXepBegin();
    CxmXepTry
    {
	rVal = Py_BuildValue("i", CxTrNheldGet(self->tr));
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
CxTreeHeldGet(CxtTreeObject *self, PyObject *args)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    int held;
    uint32_t neighbor, score, bisect, reconnectA, reconnectB;

    if (PyArg_ParseTuple(args, "i", &held) == 0)
    {
	rVal = NULL;
	goto RETURN;
    }
    if (held >= CxTrNheldGet(self->tr))
    {
	Py_INCREF(PyExc_ValueError);
	rVal = PyExc_ValueError;
	goto RETURN;
    }

    CxmXepBegin();
    CxmXepTry
    {
	CxTrHeldGet(self->tr, held, &neighbor, &score);
//	CxTrTbrNeighborGet(self->tr, neighbor,
//			    &bisect, &reconnectA, &reconnectB);

	rVal = Py_BuildValue("(i(iii))", score,
			       bisect, reconnectA, reconnectB);
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
