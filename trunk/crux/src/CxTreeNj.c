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

//#define CxmTreeNjRandomize
#define CxmTreeNjVerbose
//#define CxmTreeNjDump

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
CxpTreeNjNodesJoin(float *aD, float *aRScaled, CxtNodeObject **aNodes,
		   CxtTreeObject *aTree, long aNleft, long aXMin, long aYMin,
		   CxtNodeObject **rNode, float *rDistX, float *rDistY)
{
    CxtNodeObject *node;
    float distX, distY;
    CxtEdgeObject *edgeX, *edgeY;
    long iMin;

    /* Join the nodes that have the minimum transformed distance between
     * them. */
    node = CxNodeNew(aTree);
    edgeX = CxEdgeNew(aTree);
    CxEdgeAttach(edgeX, node, aNodes[aXMin]);
    iMin = CxpTreeNjXy2i(aNleft, aXMin, aYMin);
    distX = (aD[iMin] + aRScaled[aXMin] - aRScaled[aYMin]) / 2;
    CxEdgeLengthSet(edgeX, distX);

    edgeY = CxEdgeNew(aTree);
    CxEdgeAttach(edgeY, node, aNodes[aYMin]);
    distY = aD[iMin] - distX;
    CxEdgeLengthSet(edgeY, distY);

    *rNode = node;
    *rDistX = distX;
    *rDistY = distY;
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
	 x < aXMin;
	 x++)
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

    /* (x == aXMin) */
    iY += aNleft - 2 - x;
    x++;

    for (;
	 x < aYMin;
	 x++)
    {
	iX++;
	dist = aD[iX];
	aR[x] -= dist;
	aR[aXMin] -= dist;

	dist = aD[iY];
	iY += aNleft - 2 - x;
	aR[x] -= dist;
	aR[aYMin] -= dist;
    }

    /* (x == aYMin) */
    iX++;
    dist = aD[iX];
    aR[x] -= dist;
    aR[aXMin] -= dist;
    x++;

    for (;
	 x < aNleft;
	 x++)
    {
	iX++;
	dist = aD[iX];
	aR[x] -= dist;
	aR[aXMin] -= dist;

	iY++;
	dist = aD[iY];
	aR[x] -= dist;
	aR[aYMin] -= dist;
    }
}

CxmpInline void
CxpTreeNjCompact(float *aD, float *aR, float *aRScaled, CxtNodeObject **aNodes,
		 long aNleft, long aXMin, long aYMin, CxtNodeObject *aNode,
		 float aDistX, float aDistY)
{
    long x, iX, iY;
    float dist;

    /* Insert the new node into r. */
    aNodes[aXMin] = aNode;

    /* Calculate distances to the new node, and add them to r.  This clobbers
     * old distances, just after the last time they are needed. */
    for (x = 0,
	     iX = aXMin - 1,
	     iY = aYMin - 1;
	 x < aXMin;
	 x++)
    {
	dist = ((aD[iX] - aDistX) + (aD[iY] - aDistY)) / 2;
	aD[iX] = dist;
	iX += aNleft - 2 - x;
	iY += aNleft - 2 - x;
	aR[x] += dist;
	aR[aXMin] += dist;
    }

    /* (x == aXMin) */
    iY += aNleft - 2 - x;
    x++;

    for (;
	 x < aYMin;
	 x++)
    {
	iX++;
	dist = ((aD[iX] - aDistX) + (aD[iY] - aDistY)) / 2;
	aD[iX] = dist;
	iY += aNleft - 2 - x;
	aR[x] += dist;
	aR[aXMin] += dist;
    }

    /* (x == aYMin) */
    iX++;
    x++;

    for (;
	 x < aNleft;
	 x++)
    {
	iX++;
	iY++;
	dist = ((aD[iX] - aDistX) + (aD[iY] - aDistY)) / 2;
	aD[iX] = dist;
	aR[x] += dist;
	aR[aXMin] += dist;
    }

    /* Fill in the remaining gap (aYMin row/column), by moving the first row
     * into the gap.  The first row can be removed from the matrix in constant
     * time, whereas collapsing the gap directly would require a series of
     * memmove() calls, and leaving the gap would result in increased cache
     * misses. */
    for (x = 1,
	     iX = x - 1,
	     iY = aNleft + aYMin - 3;
	 x < aYMin;
	 x++)
    {
	aD[iY] = aD[iX];
	iY += aNleft - 2 - x;
	iX++;
    }

    /* (x == aYMin) */
    iX++;
    x++;

    for (;
	 x < aNleft;
	 x++)
    {
	iY++;
	aD[iY] = aD[iX];
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

CxmpInline bool
CxpTreeNjPairClusterOk(float *aD, float *aRScaled, long aNleft,
		       long aA, long aB)
{
    bool retval;
    long x, iA, iB;
    float distAB, dist;

    CxmAssert(aA < aB);

    /* Compare the transformed distances from {aA, aB} to every other node, and
     * make sure that the two are always closer. */
    distAB = aD[CxpTreeNjXy2i(aNleft, aA, aB)] - (aRScaled[aA] + aRScaled[aB]);

    if (aB < aNleft - 1)
    {
	for (x = aB + 1,
		 iB = CxpTreeNjXy2i(aNleft, aB, aB + 1);
	     x < aNleft;
	     x++)
	{
	    /* Make sure aA and aB are closer together than any nodes are to
	     * either aA or aB. */

	    /* Row aA was already done in CxpTreeNjCluster(). */

	    dist = aD[iB] - (aRScaled[x] + aRScaled[aB]);
	    if (dist < distAB)
	    {
		retval = false;
		goto RETURN;
	    }
	    iB++;
	}
    }

    for (x = 0,
	     iA = aA - 1,
	     iB = aB - 1;
	 x < aA;
	 x++)
    {
	/* Make sure aA and aB are closer together than any nodes are to either
	 * aA or aB. */
	dist = aD[iA] - (aRScaled[x] + aRScaled[aA]);
	if (dist < distAB)
	{
	    retval = false;
	    goto RETURN;
	}
	iA += aNleft - 2 - x;

	dist = aD[iB] - (aRScaled[x] + aRScaled[aB]);
	if (dist < distAB)
	{
	    retval = false;
	    goto RETURN;
	}
	iB += aNleft - 2 - x;
    }

    /* (x == aA) */
    iB += aNleft - 2 - x;
    x++;

    for (;
	 x < aB;
	 x++)
    {
	/* Make sure aA and aB are closer together than any nodes are to either
	 * aA or aB. */
	iA++;
	dist = aD[iA] - (aRScaled[x] + aRScaled[aA]);
	if (dist < distAB)
	{
	    retval = false;
	    goto RETURN;
	}

	dist = aD[iB] - (aRScaled[x] + aRScaled[aB]);
	if (dist < distAB)
	{
	    retval = false;
	    goto RETURN;
	}
	iB += aNleft - 2 - x;
    }

    retval = true;
#ifdef CxmTreeNjRandomize
    {
	static bool inited = false;

	if (inited == false)
	{
	    time_t t;
	    time(&t);
	    srand(t);
	    inited = true;
	}
	retval = rand() & 1;
    }
#endif
    RETURN:
    return retval;
}

static void
CxpTreeNjCluster(float **arD, float **arR, float **arRScaled,
		 CxtNodeObject ***arNodes, long *arNleft, CxtTreeObject *aTree)
{
    long x, y, min;
    float *dElm;
    float dist, minDist, distX, distY;
    float *d = *arD;
    float *r = *arR;
    float *rScaled = *arRScaled;
    CxtNodeObject *node;
    CxtNodeObject **nodes = *arNodes;
    long nleft = *arNleft;

    for (x = 0; x < nleft - 1 && nleft > 2;) /* y indexes one past x. */
    {
	/* Find the minimum distance from the node on row x to any other
	 * node that comes after this one in the matrix.  This has the effect of
	 * trying each node pairing only once. */
	if (x < nleft - 2)
	{
	    for (y = x + 1,
		     dElm = &d[CxpTreeNjXy2i(nleft, x, y)],
		     minDist = HUGE_VAL;
		 y < nleft;
		 y++)
	    {
		dist = *dElm - (rScaled[x] + rScaled[y]);
		dElm++;

		if (dist < minDist)
		{
		    minDist = dist;
		    min = y;
		}
	    }
	    CxmAssert(minDist != HUGE_VAL);
	}
	else
	{
	    min = x + 1;
	}

	if (CxpTreeNjPairClusterOk(d, rScaled, nleft, x, min))
	{
#ifdef CxmTreeNjDump
	    CxpTreeDump(d, r, rScaled, nodes, nleft);
#endif
	    CxpTreeNjNodesJoin(d, rScaled, nodes, aTree, nleft, x, min,
			       &node, &distX, &distY);
	    CxpTreeNjRSubtract(d, r, nleft, x, min);
	    CxpTreeNjCompact(d, r, rScaled, nodes, nleft, x, min,
			     node, distX, distY);
	    CxpTreeNjDiscard(&d, &r, &rScaled, &nodes, nleft);
	    nleft--;
	    CxpTreeNjRScaledUpdate(rScaled, r, nleft);

	    /* The indexing of the matrix is shifted as a result of having
             * removed the first row.  Set x such that joining this row is
             * immediately tried again.  This isn't ideal, in that this only
             * tries to join the new node with nodes that come after it.
             * However, it probably isn't worth the cache miss penalty of
             * finding the node that is closest to this one (requires iterating
             * over a column).
	     *
	     * Note that if x is 0, then the row is now at (min - 1). */
	    if (x > 0)
	    {
		x--;
	    }
	}
	else
	{
	    x++;
	}
    }

    *arD = d;
    *arR = r;
    *arRScaled = rScaled;
    *arNodes = nodes;
    *arNleft = nleft;
}

#ifdef CxmTreeNjDump
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
    long nleft;
    CxtNodeObject *node;
#ifdef CxmTreeNjVerbose
    time_t t;
    struct tm *tm;
    time_t starttime;
#endif

    CxmCheckPtr(aDistMatrix);
    CxmAssert(aNtaxa > 1);

    /* Initialize distance matrix, r, rScaled, and nodes. */
    if ((dOrig = d = CxpTreeNjMatrixInit(aDistMatrix, aNtaxa)) == NULL)
    {
	retval = true;
	goto RETURN;
    }

#ifdef CxmTreeNjVerbose
    time(&starttime);
#endif

    rOrig = r = CxpTreeNjRInit(d, aNtaxa);
    rScaledOrig = rScaled = CxpTreeNjRScaledInit(aNtaxa);
    nodesOrig = nodes = CxpTreeNjNodesInit(aTree, aNtaxa);

    nleft = aNtaxa;

    /* Iteratively try all clusterings, until only two rows are left. */
    CxpTreeNjRScaledUpdate(rScaled, r, nleft);
    while (nleft > 2)
    {
	CxpTreeNjCluster(&d, &r, &rScaled, &nodes, &nleft, aTree);
    }

    node = CxpTreeNjFinalJoin(d, nodes, aTree);

    /* Set the tree base. */
    CxTreeBaseSet(aTree, node);
#ifdef CxmTreeNjVerbose
    time(&t);
    tm = localtime(&t);
    fprintf(stderr, "%d/%02d/%02d %02d:%02d:%02d: %d second%s\n",
	    tm->tm_year + 1900, tm->tm_mon + 1, tm->tm_mday + 1,
	    tm->tm_hour, tm->tm_min, tm->tm_sec,
	    (int)(t - starttime), (t - starttime) == 1 ? "" : "s");
#endif

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
