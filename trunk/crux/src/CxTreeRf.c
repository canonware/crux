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

typedef struct
{
    unsigned nBits;
    unsigned char *bits;
} CxtTreeRfVec;

typedef struct
{
    CxtTreeRfVec treeVec;
    CxtTreeRfVec *edgeVecs;
} CxtTreeRfTree;

/******************************************************************************/

static void
CxpTreeRfVecNew(CxtTreeRfVec *aVec, unsigned aNBits)
{
    aVec->nBits = 0;
    aVec->bits = NULL;
}

static void
CxpTreeRfVecDelete(CxtTreeRfVec *aVec)
{
    if (aVec->bits != NULL)
    {
	CxmFree(aVec->bits);
    }
}

static int
CxpTreeRfVecCompare(CxtTreeRfVec *aVecA, CxtTreeRfVec *aVecB)
{
    int rVal;

    if (aVecA->nBits < aVecB->nBits)
    {
	rVal = -1;
    }
    else if (aVecA->nBits == aVecB->nBits)
    {
	rVal = memcmp(aVecA->bits, aVecB->bits, (aVecA->nBits >> 3));
    }
    else
    {
	rVal = 1;
    }

    return rVal;
}

static unsigned
CxpTreeRfVecNBitsGet(CxtTreeRfVec *aVec)
{
    return aVec->nBits;
}

static bool
CxpTreeRfVecGet(CxtTreeRfVec *aVec, uint32_t aBit, unsigned aVal)
{
    // XXX
}

static void
CxpTreeRfVecSet(CxtTreeRfVec *aVec, uint32_t aBit, bool aVal)
{
    /* Make sure that the bit vector is long enough. */
    if (aVec->nBits <= aBit)
    {
	unsigned nBytes, nBits;

	nBits = aBit + 1;
	nBytes = (nBits >> 3);
	if ((nBits & 0x7) != 0)
	{
	    nBytes++;
	}
	nBits = (nBytes >> 3);

	aVec->nBits = nBits;
	if (aVec->bits == NULL)
	{
	    aVec->bits = (unsigned char *) CxmCalloc(1, nBytes);
	}
	else
	{
	    aVec->bits = (unsigned char *) CxmRealloc(aVec->bits, nBytes);
	    /* XXX Zero-fill new bits. */
	}
    }

    /* Set bit. */
    // XXX
}

/******************************************************************************/

static void
CxpTreeRfBipartitionsInitRecurse(CxtTreeObject *self, CxtNodeObject *prevNode,
				 CxtNodeObject *curNode)
{
    CxtRingObject *firstRing, *curRing, *otherRing;
    CxtNodeObject *otherNode;
    uint32_t taxonNum;

    /* If this is a leaf node, note that this taxon exists in the tree. */
    if ((taxonNum = CxNodeTaxonNumGet(curNode)) != CxmTrNodeTaxonNone)
    {
	CxtTreeRfTree *treeAux;

	treeAux = (CxtTreeRfTree *) CxTreeAuxGet(self, CxmTreeObjectAuxRf);
	CxpTreeRfVecSet(&treeAux->treeVec, taxonNum, true);
    }

    /* Recurse. */
    firstRing = CxNodeRing(curNode);
    if (firstRing != NULL)
    {
	curRing = firstRing;
	do
	{
	    otherRing = CxRingOther(curRing);
	    otherNode = CxRingNode(otherRing);
	    if (otherNode != prevNode)
	    {
		CxpTreeRfBipartitionsInitRecurse(self, curNode, otherNode);
	    }
	    curRing = CxRingNext(curRing);
	} while (curRing != firstRing);
    }
}

static void
CxpTreeRfBipartitionsInit(CxtTreeObject *self)
{
    CxtNodeObject *base;

    base = CxTreeBaseGet(self);

    // XXX Initialize aux for self.

    // XXX

}

static float
CxpTreeRfDistanceCalc(CxtTreeObject *treeA, CxtTreeObject *treeB)
{
    float rVal;

    /* If the trees don't have precisely the same taxa, treat them as being
     * infinitely distant. */

    

    // XXX

    return rVal;
}

static void
CxpTreeRfBipartitionsCleanup(CxtTreeObject *self)
{
    // XXX Delete aux for self.
    // XXX
}

/******************************************************************************/

/* Calculate the RF distances between self and other trees, and return a tuple
 * of the corresponding distances. */
PyObject *
CxTreeRfTuple(CxtTreeObject *self, PyObject *args)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    PyObject *tuple, *other;
    unsigned i, size;
    float distance;

    if (PyArg_ParseTuple(args, "O!", &PyTuple_Type, &tuple) == 0)
    {
	rVal = NULL;
	goto RETURN;
    }

    CxmXepBegin();
    CxmXepTry
    {
	/* Make sure that all elements of the tuple are Trees. */
	size = PyTuple_Size(tuple);
	for (i = 0; i < size; i++)
	{
	    if (PyObject_IsInstance(PyTuple_GetItem(tuple, i),
				    (PyObject *) &CxtTree)
		== 0)
	    {
		CxError(CxgTreeTypeError, "Tree expected");
		goto ERROR;
	    }
	}

	/* Create a sorted list of edge-induced bipartitions for self. */
	CxpTreeRfBipartitionsInit(self);

	rVal = PyTuple_New(size);
	for (i = 0; i < size; i++)
	{
	    other = PyTuple_GetItem(tuple, i);

	    /* Create a sorted list of edge-induced bipartitions for other. */
	    CxpTreeRfBipartitionsInit((CxtTreeObject *) other);

	    /* Calculate the Robinson-Foulds distance between self and other. */
	    distance = CxpTreeRfDistanceCalc(self, (CxtTreeObject *) other);

	    /* Discard the bipartition data for other. */
	    CxpTreeRfBipartitionsCleanup((CxtTreeObject *) other);

	    /* Set the distance in rVal. */
	    PyTuple_SetItem(rVal, i, PyFloat_FromDouble(distance));
	}

	/* Discard the bipartition data for self. */
	CxpTreeRfBipartitionsCleanup(self);

	ERROR:
	rVal = NULL;
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

/* Calculate the RF distance between two trees and return the distance. */
PyObject *
CxTreeRfPair(CxtTreeObject *self, PyObject *args)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    CxtTreeObject *other;
    float distance;

    if (PyArg_ParseTuple(args, "O!", &CxtTree, (PyObject **) &other) == 0)
    {
	rVal = NULL;
	goto RETURN;
    }

    CxmXepBegin();
    CxmXepTry
    {
	/* Create a sorted list of edge-induced bipartitions for each tree. */
	CxpTreeRfBipartitionsInit(self);
	CxpTreeRfBipartitionsInit(other);

	/* Calculate the Robinson-Foulds distance between self and other. */
	distance = CxpTreeRfDistanceCalc(self, other);

	/* Discard the bipartition data. */
	CxpTreeRfBipartitionsCleanup(self);
	CxpTreeRfBipartitionsCleanup(other);

	/* Set rVal. */
	rVal = PyFloat_FromDouble(distance);
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
