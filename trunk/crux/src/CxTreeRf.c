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

// XXX This implementation assumes that taxon numbers start at 0, and are
// contiguous.

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
    unsigned nEdgeVecs;
} CxtTreeRfTree;

/******************************************************************************/

static void
CxpTreeRfVecNew(CxtTreeRfVec *aVec, unsigned aNBits)
{
    unsigned nBytes;

    nBytes = (aNBits >> 3);
    if ((aNBits & 0x7) != 0)
    {
	nBytes++;
    }
    aNBits = (nBytes << 3);

    aVec->nBits = aNBits;
    aVec->bits = (unsigned char *) CxmCalloc(1, nBytes);
}

static void
CxpTreeRfVecDelete(CxtTreeRfVec *aVec)
{
    CxmFree(aVec->bits);
}

CxmpInline int
CxpTreeRfVecCompare(const CxtTreeRfVec *aVecA, const CxtTreeRfVec *aVecB)
{
    int rVal;

    CxmAssert(aVecA->nBits == aVecB->nBits);

    rVal = memcmp(aVecA->bits, aVecB->bits, (aVecA->nBits >> 3));

    return rVal;
}

/* qsort(3)-compatible vector comparison function. */
static int
CxpTreeRfVecQsortCompare(const void *aA, const void *aB)
{
    CxtTreeRfVec *aVecA = (CxtTreeRfVec *) aA;
    CxtTreeRfVec *aVecB = (CxtTreeRfVec *) aB;

    return CxpTreeRfVecCompare(aVecA, aVecB);
}

static bool
CxpTreeRfVecGet(CxtTreeRfVec *aVec, uint32_t aBit)
{
    bool rVal;
    unsigned byteOffset, bitOffset;
    unsigned char byte;

    CxmAssert(aBit < aVec->nBits);

    byteOffset = aBit >> 3;
    bitOffset = aBit & 0x7;

    byte = aVec->bits[byteOffset];
    byte >>= bitOffset;
    byte &= 0x1;
    rVal = byte;

    return rVal;
}

static void
CxpTreeRfVecSet(CxtTreeRfVec *aVec, uint32_t aBit, bool aVal)
{
    unsigned byteOffset, bitOffset;
    unsigned char byte, bit, mask;

    CxmAssert(aBit < aVec->nBits);

    byteOffset = aBit >> 3;
    bitOffset = aBit & 0x7;

    byte = aVec->bits[byteOffset];

    mask = (0xff ^ (0x1 << bitOffset));
    byte &= mask;
    bit = (aVal << bitOffset);
    byte |= bit;

    aVec->bits[byteOffset] = byte;
}

static void
CxpTreeRfVecInvert(CxtTreeRfVec *aVec)
{
    unsigned i, nBytes;

    for (i = 0, nBytes = (aVec->nBits >> 3); i < nBytes; i++)
    {
	aVec->bits[i] = (~aVec->bits[i]);
    }
}

/* Create the union of aA and aB, and store the result in aA. */
static void
CxpTreeRfVecUnion(CxtTreeRfVec *aA, CxtTreeRfVec *aB)
{
    unsigned i, nBytes;

    CxmAssert(aA->nBits == aB->nBits);

    for (i = 0, nBytes = (aA->nBits >> 3); i < nBytes; i++)
    {
	aA->bits[i] |= aB->bits[i];
    }
}

/******************************************************************************/

static CxtTreeRfVec *
CxpTreeRfBipartitionsInitRecurse(CxtTreeObject *self,
				 CxtTreeRfTree *aRfTree,
				 CxtNodeObject *prevNode,
				 CxtNodeObject *curNode)
{
    CxtTreeRfVec *rVal, *vec;
    CxtRingObject *firstRing, *curRing, *otherRing;
    CxtNodeObject *otherNode;
    unsigned ntaxa;
    uint32_t taxonNum;

    ntaxa = CxTreeNtaxaGet(self);
    rVal = &aRfTree->edgeVecs[aRfTree->nEdgeVecs];
    aRfTree->nEdgeVecs++;
    CxpTreeRfVecNew(rVal, ntaxa);

    if ((taxonNum = CxNodeTaxonNumGet(curNode)) != CxmTrNodeTaxonNone)
    {
	/* Leaf node. */

	/* Note that this taxon exists in the tree. */
	CxpTreeRfVecSet(&aRfTree->treeVec, taxonNum, true);

	/* Mark this taxon as being in the set. */
	CxpTreeRfVecSet(rVal, taxonNum, true);
    }
    else
    {
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
		    vec = CxpTreeRfBipartitionsInitRecurse(self, aRfTree,
							   curNode, otherNode);
		    CxpTreeRfVecUnion(rVal, vec);
		}
		curRing = CxRingNext(curRing);
	    } while (curRing != firstRing);
	}
    }

    return rVal;
}

static CxtTreeRfTree *
CxpTreeRfBipartitionsInit(CxtTreeObject *self)
{
    CxtTreeRfTree *rVal;
    CxtNodeObject *baseNode, *otherNode;
    CxtRingObject *firstRing, *curRing, *otherRing;
    unsigned ntaxa, i;
    uint32_t taxonNum;

    ntaxa = CxTreeNtaxaGet(self);

    /* Initialize rVal. */
    rVal = (CxtTreeRfTree *) CxmMalloc(sizeof(CxtTreeRfTree));
    CxpTreeRfVecNew(&rVal->treeVec, ntaxa);
    /* Allocate at least enough slots in edgeVecs. */
    rVal->edgeVecs = (CxtTreeRfVec *) CxmCalloc(ntaxa * 2,
						sizeof(CxtTreeRfVec));
    rVal->nEdgeVecs = 0;

    /* Recurse through tree and determine bipartitions. */
    baseNode = CxTreeBaseGet(self);
    if (baseNode != NULL)
    {
	if ((taxonNum = CxNodeTaxonNumGet(baseNode)) != CxmTrNodeTaxonNone)
	{
	    /* Leaf node. */

	    /* Note that this taxon exists in the tree. */
	    CxpTreeRfVecSet(&rVal->treeVec, taxonNum, true);
	}


	firstRing = CxNodeRing(baseNode);
	if (firstRing != NULL)
	{
	    curRing = firstRing;
	    do
	    {
		otherRing = CxRingOther(curRing);
		otherNode = CxRingNode(otherRing);
		CxpTreeRfBipartitionsInitRecurse(self, rVal,
						 baseNode, otherNode);
		curRing = CxRingNext(curRing);
	    } while (curRing != firstRing);
	}
    }

    /* Iterate through bipartitions and make sure that taxon 0 is never in the
     * set represented by the bit vector. */
    for (i = 0; i < rVal->nEdgeVecs; i++)
    {
	if (CxpTreeRfVecGet(&rVal->edgeVecs[i], 0))
	{
	    CxpTreeRfVecInvert(&rVal->edgeVecs[i]);
	}
    }

    /* Sort the list of bipartions. */
    qsort(rVal->edgeVecs, rVal->nEdgeVecs, sizeof(CxtTreeRfVec),
	  CxpTreeRfVecQsortCompare);

    return rVal;
}

static float
CxpTreeRfDistanceCalc(CxtTreeObject *aTreeA, CxtTreeRfTree *aTreeARfVec,
		      CxtTreeObject *aTreeB, CxtTreeRfTree *aTreeBRfVec)
{
    float rVal;
    unsigned iA, iB, nUniqueA, nUniqueB;
    int relation;

    /* If the trees don't have precisely the same taxa, treat them as being
     * infinitely distant. */
    if (CxpTreeRfVecCompare(&aTreeARfVec->treeVec, &aTreeBRfVec->treeVec) != 0)
    {
	rVal = 1.0;
	goto RETURN;
    }

    /* Count the number of unique bipartitions in each tree. */
    for (iA = iB = nUniqueA = nUniqueB = 0;
	 iA < aTreeARfVec->nEdgeVecs && iB < aTreeBRfVec->nEdgeVecs;
	 )
    {
	relation = CxpTreeRfVecCompare(&aTreeARfVec->edgeVecs[iA],
				       &aTreeBRfVec->edgeVecs[iB]);
	switch (relation)
	{
	    case -1:
	    {
		/* aTreeA has a unique bipartion. */
		nUniqueA++;
		iA++;
		break;
	    }
	    case 0:
	    {
		iA++;
		iB++;
		break;
	    }
	    case 1:
	    {
		/* aTreeB has a unique bipartion. */
		nUniqueB++;
		iB++;
		break;
	    }
	    default:
	    {
		CxmNotReached();
	    }
	}
    }

    if (iA < aTreeARfVec->nEdgeVecs)
    {
	nUniqueA += aTreeARfVec->nEdgeVecs - iA;
    }
    else if (iB < aTreeBRfVec->nEdgeVecs)
    {
	nUniqueB += aTreeBRfVec->nEdgeVecs - iB;
    }

    /* Convert counts to the Robinson-Foulds distance. */
    rVal = ((((float) nUniqueA / (float) aTreeARfVec->nEdgeVecs)
	     + ((float) nUniqueB / (float) aTreeBRfVec->nEdgeVecs))
	    / 2.0);

    RETURN:
    return rVal;
}

static void
CxpTreeRfBipartitionsCleanup(CxtTreeObject *self, CxtTreeRfTree *aRfTree)
{
    unsigned i;

    CxpTreeRfVecDelete(&aRfTree->treeVec);

    for (i = 0; i < aRfTree->nEdgeVecs; i++)
    {
	CxpTreeRfVecDelete(&aRfTree->edgeVecs[i]);
    }
    CxmFree(aRfTree->edgeVecs);

    CxmFree(aRfTree);
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
    CxtTreeRfTree *rfTreeSelf, *rfTreeOther;
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
	rfTreeSelf = CxpTreeRfBipartitionsInit(self);

	rVal = PyTuple_New(size);
	for (i = 0; i < size; i++)
	{
	    other = PyTuple_GetItem(tuple, i);

	    if ((CxtTreeObject *) other == self)
	    {
		distance = 0.0;
	    }
	    else
	    {
		/* Create a sorted list of edge-induced bipartitions for
		 * other. */
		rfTreeOther =
		    CxpTreeRfBipartitionsInit((CxtTreeObject *) other);

		/* Calculate the Robinson-Foulds distance between self and
		 * other. */
		distance =
		    CxpTreeRfDistanceCalc(self, rfTreeSelf,
					  (CxtTreeObject *) other,
					  rfTreeOther);

		/* Discard the bipartition data for other. */
		CxpTreeRfBipartitionsCleanup((CxtTreeObject *) other,
					     rfTreeOther);
	    }

	    /* Set the distance in rVal. */
	    PyTuple_SetItem(rVal, i, PyFloat_FromDouble(distance));
	}

	/* Discard the bipartition data for self. */
	CxpTreeRfBipartitionsCleanup(self, rfTreeSelf);

	break;
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
    CxtTreeRfTree *rfTreeSelf, *rfTreeOther;
    float distance;

    if (PyArg_ParseTuple(args, "O!", &CxtTree, (PyObject **) &other) == 0)
    {
	rVal = NULL;
	goto RETURN;
    }
    if (other == self)
    {
	rVal = PyFloat_FromDouble(0.0);
	goto RETURN;
    }

    CxmXepBegin();
    CxmXepTry
    {
	/* Create a sorted list of edge-induced bipartitions for each tree. */
	rfTreeSelf = CxpTreeRfBipartitionsInit(self);
	rfTreeOther = CxpTreeRfBipartitionsInit(other);

	/* Calculate the Robinson-Foulds distance between self and other. */
	distance = CxpTreeRfDistanceCalc(self, rfTreeSelf,
					 other, rfTreeOther);

	/* Discard the bipartition data. */
	CxpTreeRfBipartitionsCleanup(self, rfTreeSelf);
	CxpTreeRfBipartitionsCleanup(other, rfTreeOther);

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
