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

typedef struct CxsTreeNjr CxtTreeNjr;

/* The following function can be used to convert from row/column matrix
 * coordinates to array offsets for neighbor-joining:
 *
 *   n : Number of nodes currently in the matrix.
 *   x : Row.
 *   y : Column.
 *
 *                        2
 *                       x  + 3x
 *   f(n,x,y) = nx + y - ------- - 1
 *                          2
 */
CxmpInline unsigned long
CxpTreeNjXy2i(unsigned long aN, unsigned long aX, unsigned long aY)
{
    CxmAssert(aX < aN);
    CxmAssert(aY < aN);
    CxmAssert(aX < aY);

    return aN * aX + aY - (((aX + 3) * aX) >> 1) - 1;
}

static float *
CxpTreeNjMatrixInit(PyObject *aDistMatrix, long aNtaxa)
{
    float *retval;
    float *dElm;
    long x, y;
    PyObject *result;

    /* Allocate an array that is large enough to hold the distances. */
    retval = (float *) CxmMalloc(sizeof(float)
				 * (CxpTreeNjXy2i(aNtaxa, aNtaxa - 2,
						  aNtaxa - 1)
				    + 1));

    /* Initialize untransformed distances. */
    for (x = 0, dElm = retval; x < aNtaxa; x++)
    {
	for (y = x + 1; y < aNtaxa; y++)
	{
	    result = PyEval_CallMethod(aDistMatrix, "distanceGet", "(ll)",
				       x, y);
	    if (PyFloat_Check(result))
	    {
		*dElm = (float) PyFloat_AsDouble(result);
		Py_DECREF(result);
	    }
	    else if (PyInt_Check(result))
	    {
		*dElm = (float) PyInt_AsLong(result);
		Py_DECREF(result);
	    }
	    else
	    {
		Py_DECREF(result);
		CxError(CxgTreeTypeError,
			"Int or float distance expected (%ld, %ld)",
			x, y);
		CxmFree(retval);
		retval = NULL;
		goto RETURN;
	    }

	    dElm++;
	}
    }

    RETURN:
    return retval;
}

static float *
CxpTreeNjRInit(float *aD, long aNtaxa)
{
    float *retval;
    float dist, *dElm;
    long x, y;

    retval = (float *) CxmMalloc(sizeof(float) * aNtaxa);
    
    /* Calculate r (sum of distances to other nodes) for each node. */
    for (x = 0; x < aNtaxa; x++)
    {
	retval[x] = 0.0;
    }

    for (x = 0, dElm = aD; x < aNtaxa; x++)
    {
	for (y = x + 1; y < aNtaxa; y++)
	{
	    dist = *dElm;
	    dElm++;

	    retval[x] += dist;
	    retval[y] += dist;
	}
    }

    return retval;
}

static float *
CxpTreeNjRScaledInit(long aNtaxa)
{
    float *retval;

    retval = (float *) CxmMalloc(sizeof(float) * aNtaxa);

    return retval;
}

static CxtNodeObject **
CxpTreeNjNodesInit(CxtTreeObject *aTree, long aNtaxa)
{
    CxtNodeObject **retval;
    long x;

    retval = (CxtNodeObject **) CxmMalloc(sizeof(CxtNodeObject *) * aNtaxa);

    /* Create a node for each taxon in the matrix. */
    for (x = 0; x < aNtaxa; x++)
    {
	retval[x] = CxNodeNew(aTree);
	CxNodeTaxonNumSet(retval[x], x);
    }

    return retval;
}

CxmpInline void
CxpTreeNjRScaledUpdate(float *aRScaled, float *aR, long aNleft)
{
    long x;

    /* Calculate rScaled (r/(nleft-2)) for each node. */
    for (x = 0; x < aNleft; x++)
    {
	aRScaled[x] = aR[x] / (aNleft - 2);
    }
}

CxmpInline void
CxpTreeNjMinFind(float *aD, float *aRScaled, long aNleft,
		 long *rXMin, long *rYMin, long *rIMin)
{
    float *dElm, transMin, transCur;
    long x, y, xMin, yMin, iMin;

    /* Calculate the transformed distance for each pairwise distance.  Keep
     * track of the minimum transformed distance, so that the corresponding
     * nodes can be joined.  Ties are broken arbitrarily (the first minimum
     * found is used).
     *
     * This is by far the most time-consuming portion of neighbor joining.
     *
     * Use pointer arithmetic (dElm), rather than d[i].  This appears to reduce
     * register pressure on x86, and has a significant positive performance
     * impact. */
#ifdef CxmCcSilence
    xMin = yMin = 0;
#endif
    for (x = 0, dElm = aD, transMin = HUGE_VAL; x < aNleft; x++)
    {
	for (y = x + 1; y < aNleft; y++)
	{
	    transCur = *dElm - (aRScaled[x] + aRScaled[y]);
	    dElm++;

	    if (transCur < transMin)
	    {
		xMin = x;
		yMin = y;
		transMin = transCur;
	    }
	}
    }
    iMin = CxpTreeNjXy2i(aNleft, xMin, yMin);

    *rXMin = xMin;
    *rYMin = yMin;
    *rIMin = iMin;
}

CxmpInline CxtNodeObject *
CxpTreeNjNodesJoin(float *aD, float *aRScaled, CxtNodeObject **aNodes,
		   CxtTreeObject *aTree, long aXMin, long aYMin, long aIMin)
{
    CxtNodeObject *retval;
    float distX, distY;
    CxtEdgeObject *edgeX, *edgeY;
    
    /* Join the nodes that have the minimum transformed distance between
     * them. */
    retval = CxNodeNew(aTree);
    edgeX = CxEdgeNew(aTree);
    CxEdgeAttach(edgeX, retval, aNodes[aXMin]);
    distX = (aD[aIMin] + aRScaled[aXMin] - aRScaled[aYMin]) / 2;
    CxEdgeLengthSet(edgeX, distX);

    edgeY = CxEdgeNew(aTree);
    CxEdgeAttach(edgeY, retval, aNodes[aYMin]);
    distY = aD[aIMin] - distX;
    CxEdgeLengthSet(edgeY, distY);

    return retval;
}

CxmpInline void
CxpTreeNjRSubtract(float *aD, float *aR, long aNleft, long aXMin, long aYMin)
{
    long x, iX, iY;
    float dist;

    /* Subtract old distances from r. */
    for (x = 0,
	     iX = aXMin - 1,
	     iY = aYMin - 1;
	 x < aNleft;
	 x++)
    {
	if (x < aXMin)
	{
	    dist = aD[iX];
	    iX += aNleft - 2 - x;
	    aR[x] -= dist;
	    aR[aXMin] -= dist;

	    dist = aD[iY];
	    iY += aNleft - 2 - x;
	    aR[x] -= dist;
	    aR[aYMin] -= dist;
	}
	else if (x > aXMin)
	{
	    iX++;
	    dist = aD[iX];
	    aR[x] -= dist;
	    aR[aXMin] -= dist;

	    if (x < aYMin)
	    {
		dist = aD[iY];
		iY += aNleft - 2 - x;
		aR[x] -= dist;
		aR[aYMin] -= dist;
	    }
	    else if (x > aYMin)
	    {
		iY++;
		dist = aD[iY];
		aR[x] -= dist;
		aR[aYMin] -= dist;
	    }
	}
	else // (x == aXMin)
	{
	    iY += aNleft - 2 - x;
	}
    }
}

CxmpInline void
CxpTreeNjCompact(float *aD, float *aR, float *aRScaled, CxtNodeObject **aNodes,
		 long aNleft, long aXMin, long aYMin, CxtNodeObject *aNode)
{
    long x, iX, iY;
    float dist, distX, distY;
#ifdef CxmCcSilence
    distX = distY = 0.0;
#endif

    /* Insert the new node into r. */
    aNodes[aXMin] = aNode;

    /* Calculate distances to the new node, and add them to r.  This clobbers
     * old distances, just after the last time they are needed. */
    for (x = 0,
	     iX = aXMin - 1,
	     iY = aYMin - 1;
	 x < aNleft;
	 x++)
    {
	if (x < aXMin)
	{
	    dist = ((aD[iX] - distX) + (aD[iY] - distY)) / 2;
	    aD[iX] = dist;
	    iX += aNleft - 2 - x;
	    iY += aNleft - 2 - x;
	    aR[x] += dist;
	    aR[aXMin] += dist;
	}
	else if (x > aXMin)
	{
	    if (x < aYMin)
	    {
		iX++;
		dist = ((aD[iX] - distX) + (aD[iY] - distY)) / 2;
		aD[iX] = dist;
		iY += aNleft - 2 - x;
		aR[x] += dist;
		aR[aXMin] += dist;
	    }
	    else if (x > aYMin)
	    {
		iX++;
		iY++;
		dist = ((aD[iX] - distX) + (aD[iY] - distY)) / 2;
		aD[iX] = dist;
		aR[x] += dist;
		aR[aXMin] += dist;
	    }
	    else // if (x == aYMin)
	    {
		iX++;
	    }
	}
	else // if (x == aXMin)
	{
	    iY += aNleft - 2 - x;
	}
    }

    /* Fill in the remaining gap (aYMin row/column), by moving the first row
     * into the gap.  The first row can be removed from the matrix in constant
     * time, whereas collapsing the gap directly would require a series of
     * memmove() calls, and leaving the gap would result in increased cache
     * misses. */
    for (x = 1,
	     iX = x - 1,
	     iY = aNleft + aYMin - 3;
	 x < aNleft;
	 x++)
    {
	if (x < aYMin)
	{
	    aD[iY] = aD[iX];
	    iY += aNleft - 2 - x;
	}
	else if (x > aYMin)
	{
	    iY++;
	    aD[iY] = aD[iX];
	}
	iX++;
    }
    /* Fill in the gap in r, aRScaled, and nodes. */
    aR[aYMin] = aR[0];
    aRScaled[aYMin] = aRScaled[0];
    aNodes[aYMin] = aNodes[0];
}

CxmpInline void
CxpTreeNjDiscard(float **arD, float **arR, float **arRScaled,
		 CxtNodeObject ***arNodes, long aNleft)
{
    /* Move pointers forward, which removes the first row. */
    *arD = &(*arD)[aNleft - 1];
    *arR = &(*arR)[1];
    *arRScaled = &(*arRScaled)[1];
    *arNodes = &(*arNodes)[1];
}

static CxtNodeObject *
CxpTreeNjFinalJoin(float *aD, CxtNodeObject **aNodes, CxtTreeObject *aTree)
{
    CxtEdgeObject *edge;
    
    /* Join the remaining two nodes. */
    edge = CxEdgeNew(aTree);
    CxEdgeAttach(edge, aNodes[0], aNodes[1]);
    CxEdgeLengthSet(edge, aD[0]);

    return aNodes[0];
}

//#define CxmTreeNjVerbose
#ifdef CxmTreeNjVerbose
static void
CxpTreeDump(float *aD, float *aR, float *aRScaled, CxtNodeObject **aNodes,
	    long aNleft)
{
    PyObject *result;
    float *dElm;
    long x, y;
     
    fprintf(stderr,
	    "----------------------------------------"
	    "----------------------------------------\n");
    for (x = 0, dElm = aD; x < aNleft; x++)
    {
	fprintf(stderr, "%*s", (int) x * 9 + (!!x), " ");
	for (y = x + 1; y < aNleft; y++)
	{
	    fprintf(stderr, " %8.4f", *dElm);
	    dElm++;
	}
	result = CxNodeTaxonNumGet(aNodes[x]);
	if (result != Py_None)
	{
	    fprintf(stderr, " || %8.4f %8.4f (node %ld %p)\n",
		    aR[x], aRScaled[x],
		    PyInt_AsLong(result), aNodes[x]);
	}
	else
	{
	    fprintf(stderr, " || %8.4f %8.4f (internal node %p)\n",
		    aR[x], aRScaled[x], aNodes[x]);
	}
    }
}
#endif

/*
 * Create a tree from a pairwise distance matrix, using the neighbor-joining
 * algorithm.
 *
 * The matrix is actually stored as an upper-triangle symmetric matrix, with
 * additional bookkeeping, as necessary for the neighbor-joining algorithm.
 * For example (m is current matrix size):
 *
 *                             |  r  ||
 *                             | --- ||
 *   | A | B | C | D | E ||  r | n-2 ||
 *   +===+===+===+===+===++====+=====++===
 *       | 0 | 1 | 2 | 3 ||  6 | 2.0 || A
 *       +---+---+---+---++----+-----++---
 *           | 4 | 5 | 6 || 15 | 5.0 || B
 *           +---+---+---++----+-----++---
 *               | 7 | 8 || 20 | 6.7 || C
 *               +---+---++----+-----++---
 *                   | 9 || 23 | 7.7 || D
 *                   +---++----+-----++---
 *                       || 26 | 8.7 || E
 *                       ++----+-----++---
 *
 * is stored as:
 *
 *   d:
 *   /---+---+---+---+---+---+---+---+---+---\
 *   | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
 *   \---+---+---+---+---+---+---+---+---+---/
 *
 *   r:
 *   /------+------+------+------+------\
 *   |  6   | 15   | 20   | 23   | 26   |
 *   \------+------+------+------+------/
 *
 *   rScaled:
 *   /------+------+------+------+------\
 *   |  2.0 |  5.0 |  6.7 |  7.7 |  8.7 |
 *   \------+------+------+------+------/
 *
 *   nodes:
 *   /------+------+------+------+------\
 *   |  A   |  B   |  C   |  D   |  E   |
 *   \------+------+------+------+------/
 */
static bool
CxpTreeNj(CxtTreeObject *aTree, PyObject *aDistMatrix, long aNtaxa)
{
    bool retval;
    float *dOrig, *d; /* Distance matrix. */
    float *rOrig, *r; /* Distance sums. */
    float *rScaledOrig, *rScaled; /* Scaled distance sums: r/(nleft-2)). */
    CxtNodeObject **nodesOrig, **nodes; /* Nodes associated with each row. */
    long nleft, iMin, xMin, yMin;
    CxtNodeObject *node;

    CxmCheckPtr(aDistMatrix);
    CxmAssert(aNtaxa > 1);

    /* Initialize distance matrix, r, rScaled, and nodes. */
    if ((dOrig = d = CxpTreeNjMatrixInit(aDistMatrix, aNtaxa)) == NULL)
    {
	retval = true;
	goto RETURN;
    }
    rOrig = r = CxpTreeNjRInit(d, aNtaxa);
    rScaledOrig = rScaled = CxpTreeNjRScaledInit(aNtaxa);
    nodesOrig = nodes = CxpTreeNjNodesInit(aTree, aNtaxa);

    /* Iteratitively join two nodes in the matrix, until only two are left. */
    for (nleft = aNtaxa; nleft > 2; nleft--)
    {
#if (0)
	{
	    static bool initialized = false;
	    time_t starttime, t;
	    long ntaxa;
	    struct tm *tm;

	    if (initialized == false)
	    {
		ntaxa = nleft;
		time(&starttime);
		initialized = true;
	    }

	    time(&t);
	    tm = localtime(&t);
	    if (nleft < ntaxa)
	    {
		fprintf(stderr,
			"%d/%02d/%02d %02d:%02d:%02d [%ld]: %1.3f sec/join\n",
			tm->tm_year + 1900, tm->tm_mon + 1, tm->tm_mday + 1,
			tm->tm_hour, tm->tm_min, tm->tm_sec, nleft,
			((float) (t - starttime)) / ((float) (ntaxa - nleft)));
	    }
	}
#endif
#ifdef CxmTreeNjVerbose
	CxpTreeDump(d, r, rScaled, nodes, nleft);
#endif

	CxpTreeNjRScaledUpdate(rScaled, r, nleft);
	CxpTreeNjMinFind(d, rScaled, nleft, &xMin, &yMin, &iMin);
	node = CxpTreeNjNodesJoin(d, rScaled, nodes, aTree, xMin, yMin, iMin);
	CxpTreeNjRSubtract(d, r, nleft, xMin, yMin);
	CxpTreeNjCompact(d, r, rScaled, nodes, nleft, xMin, yMin, node);
	CxpTreeNjDiscard(&d, &r, &rScaled, &nodes, nleft);
    }

    node = CxpTreeNjFinalJoin(d, nodes, aTree);

    /* Set the tree base. */
    CxTreeBaseSet(aTree, node);

    retval = false;
    /* Clean up. */
    CxmFree(nodesOrig);
    CxmFree(rScaledOrig);
    CxmFree(rOrig);
    CxmFree(dOrig);
    RETURN:
    return retval;
}

PyObject *
CxTreeNj(CxtTreeObject *self, PyObject *args)
{
    PyObject *retval, *result, *distMatrix;
    long ntaxa;

    if (PyArg_ParseTuple(args, "O", &distMatrix) == 0)
    {
	retval = NULL;
	goto RETURN;
    }

    result = PyEval_CallMethod(distMatrix, "ntaxaGet", "()");
    if (PyInt_Check(result) == false)
    {
	CxError(CxgDistMatrixTypeError,
		"Integer expected from distMatrix.ntaxaGet()");
	retval = NULL;
	goto RETURN;
    }
    ntaxa = PyInt_AsLong(result);
    Py_DECREF(result);

    Py_INCREF(Py_None);
    retval = Py_None;
    CxmXepBegin();
    CxmXepTry
    {
	CxtTrNode oldTrNode, trNode;
	CxtNodeObject *node;

	oldTrNode = CxTrBaseGet(self->tr);

	/* Neighbor-join. */
	if (CxpTreeNj(self, distMatrix, ntaxa))
	{
	    /* Error during neighbor join. */
	    Py_DECREF(retval);
	    retval = NULL;
	    break;
	}

	/* Reference new base. */
	trNode = CxTrBaseGet(self->tr);
	if (trNode != CxmTrNodeNone)
	{
	    node = (CxtNodeObject *) CxTrNodeAuxGet(self->tr, trNode);
	    Py_INCREF(node);
	}

	/* Decref old base. */
	if (oldTrNode != CxmTrNodeNone)
	{
	    node = (CxtNodeObject *) CxTrNodeAuxGet(self->tr, oldTrNode);
	    Py_DECREF(node);
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
