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
//#define CxmTreeNjVerbose
static bool
CxpTreeNj(CxtTreeObject *aTree, PyObject *aDistMatrix, long aNtaxa)
{
    bool retval;
    float *dOrig, *d; /* Distance matrix. */
    float dist, *dElm;
    float *rOrig, *r; /* Distance sums. */
    float *rScaledOrig, *rScaled; /* Scaled distance sums: r/(nleft-2)). */
    CxtNodeObject **nodesOrig, **nodes; /* Nodes associated with each row. */
    long nleft, x, y, iX, iY, iMin, xMin, yMin;
    float transMin, transCur;
    float distX, distY;
    CxtNodeObject *node;
    CxtEdgeObject *edgeX, *edgeY, *edge;
    PyObject *result;

    CxmCheckPtr(aDistMatrix);
    CxmAssert(aNtaxa > 1);

    /* Allocate an array that is large enough to hold the distances. */
    dOrig = d = (float *) CxmMalloc(sizeof(float)
				    * (CxpTreeNjXy2i(aNtaxa, aNtaxa - 2,
						     aNtaxa - 1)
				       + 1));

    /* Allocate arrays that are large enough to hold all the distance sums and
     * nodes. */
    rOrig = r = (float *) CxmMalloc(sizeof(float) * aNtaxa);
    rScaledOrig = rScaled = (float *) CxmMalloc(sizeof(float) * aNtaxa);
    nodesOrig = nodes = (CxtNodeObject **) CxmMalloc(sizeof(CxtNodeObject)
						     * aNtaxa);

    /* Initialize untransformed distances. */
    for (x = 0, dElm = d; x < aNtaxa; x++)
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
		retval = true;
		goto RETURN;
	    }

	    dElm++;
	}
    }

    /* Create a node for each taxon in the matrix. */
    for (x = 0; x < aNtaxa; x++)
    {
	nodes[x] = CxNodeNew(aTree);
	CxNodeTaxonNumSet(nodes[x], x);
    }

    /* Calculate r (sum of distances to other nodes) for each node. */
    for (x = 0; x < aNtaxa; x++)
    {
	r[x] = 0.0;
    }

    for (x = 0, dElm = d; x < aNtaxa; x++)
    {
	for (y = x + 1; y < aNtaxa; y++)
	{
	    dist = *dElm;
	    dElm++;

	    r[x] += dist;
	    r[y] += dist;
	}
    }

    /* Iteratitively join two nodes in the matrix, until only two are left. */
    for (nleft = aNtaxa; nleft > 2; nleft--)
    {
#ifdef CxmTreeNjVerbose
	{
	    time_t t;
	    struct tm *tm;

	    time(&t);
	    tm = localtime(&t);
	    fprintf(stderr, "%d/%02d/%02d %02d:%02d:%02d [%ld]\n",
		    tm->tm_year + 1900, tm->tm_mon + 1, tm->tm_mday + 1,
		    tm->tm_hour, tm->tm_min, tm->tm_sec, nleft);
	}
	fprintf(stderr, "Size: %ld\n", nleft);
#endif

	/* Calculate rScaled (r/(nleft-2)) for each node. */
	for (x = 0; x < nleft; x++)
	{
	    rScaled[x] = r[x] / (nleft - 2);
	}

	/* Calculate the transformed distance for each pairwise distance.  Keep
	 * track of the minimum transformed distance, so that the corresponding
	 * nodes can be joined.  Ties are broken arbitrarily (the first minimum
	 * found is used).
	 *
	 * This is by far the most time-consuming portion of neighbor joining.
	 *
	 * Use pointer arithmetic (dElm), rather than d[i].  This appears to
	 * reduce register pressure on x86, and has a significant positive
	 * performance impact. */
#ifdef CxmCcSilence
	xMin = yMin = 0;
#endif
	for (x = 0, dElm = d, transMin = HUGE_VAL; x < nleft; x++)
	{
	    for (y = x + 1; y < nleft; y++)
	    {
		transCur = *dElm - (rScaled[x] + rScaled[y]);
		dElm++;

		if (transCur < transMin)
		{
		    xMin = x;
		    yMin = y;
		    transMin = transCur;
		}
	    }
	}
	iMin = CxpTreeNjXy2i(nleft, xMin, yMin);

#ifdef CxmTreeNjVerbose
	{
	    PyObject *result;
     
	    fprintf(stderr, "Min: %8.4f at d[%ld] (%ld, %ld)\n",
		    transMin, iMin, xMin, yMin);
	    fprintf(stderr,
		    "----------------------------------------"
		    "----------------------------------------\n");
	    for (x = 0, dElm = d; x < nleft; x++)
	    {
		fprintf(stderr, "%*s", (int) x * 9 + (!!x), " ");
		for (y = x + 1; y < nleft; y++)
		{
		    fprintf(stderr, " %8.4f", *dElm);
		    dElm++;
		}
		result = CxNodeTaxonNumGet(nodes[x]);
		if (result != Py_None)
		{
		    fprintf(stderr, " || %8.4f %8.4f (node %ld %p)\n",
			    r[x], rScaled[x],
			    PyInt_AsLong(result), nodes[x]);
		}
		else
		{
		    fprintf(stderr, " || %8.4f %8.4f (internal node %p)\n",
			    r[x], rScaled[x], nodes[x]);
		}
	    }
	}
#endif

	/* Join the nodes that have the minimum transformed distance between
	 * them. */
	node = CxNodeNew(aTree);
#ifdef CxmTreeNjVerbose
	fprintf(stderr, "New node %p\n", node);
#endif
	edgeX = CxEdgeNew(aTree);
	CxEdgeAttach(edgeX, node, nodes[xMin]);
	distX = (d[iMin] + rScaled[xMin] - rScaled[yMin]) / 2;
	CxEdgeLengthSet(edgeX, distX);
#ifdef CxmTreeNjVerbose
	{
	    PyObject *result = CxNodeTaxonNumGet(nodes[xMin]);

	    if (result != Py_None)
	    {
		fprintf(stderr, "  Join node %ld %p (len %.4f)\n",
			PyInt_AsLong(result), nodes[xMin],
			distX);
	    }
	    else
	    {
		fprintf(stderr, "  Join internal node %p (len %.4f)\n",
			nodes[xMin], distX);
	    }
	}
#endif

	edgeY = CxEdgeNew(aTree);
	CxEdgeAttach(edgeY, node, nodes[yMin]);
	distY = d[iMin] - distX;
	CxEdgeLengthSet(edgeY, distY);
#ifdef CxmTreeNjVerbose
	{
	    PyObject *result = CxNodeTaxonNumGet(nodes[yMin]);

	    if (result != Py_None)
	    {
		fprintf(stderr, "  Join node %ld %p (len %.4f)\n",
			PyInt_AsLong(result), nodes[yMin],
			distY);
	    }
	    else
	    {
		fprintf(stderr, "  Join internal node %p (len %.4f)\n",
			nodes[yMin], distY);
	    }
	}
#endif

	/* Subtract old distances from r. */
	for (x = 0,
		 iX = CxpTreeNjXy2i(nleft, x, xMin),
		 iY = CxpTreeNjXy2i(nleft, x, yMin);
	     x < nleft;
	     x++)
	{
	    if (x < xMin)
	    {
		dist = d[iX];
		iX += nleft - 2 - x;
		r[x] -= dist;
		r[xMin] -= dist;

		dist = d[iY];
		iY += nleft - 2 - x;
		r[x] -= dist;
		r[yMin] -= dist;
	    }
	    else if (x > xMin)
	    {
		iX++;
		dist = d[iX];
		r[x] -= dist;
		r[xMin] -= dist;

		if (x < yMin)
		{
		    dist = d[iY];
		    iY += nleft - 2 - x;
		    r[x] -= dist;
		    r[yMin] -= dist;
		}
		else if (x > yMin)
		{
		    iY++;
		    dist = d[iY];
		    r[x] -= dist;
		    r[yMin] -= dist;
		}
	    }
	    else // (x == xMin)
	    {
		iY += nleft - 2 - x;
	    }
	}

	/* Compact matrix. */

	/* Insert the new node into r. */
	nodes[xMin] = node;

	/* Calculate distances to the new node, and add them to r.  This
	 * clobbers old distances, just after the last time they are needed. */
	for (x = 0,
		 iX = CxpTreeNjXy2i(nleft, x, xMin),
		 iY = CxpTreeNjXy2i(nleft, x, yMin);
	     x < nleft;
	     x++)
	{
	    if (x < xMin)
	    {
		dist = ((d[iX] - distX) + (d[iY] - distY)) / 2;
		d[iX] = dist;
		iX += nleft - 2 - x;
		iY += nleft - 2 - x;
		r[x] += dist;
		r[xMin] += dist;
	    }
	    else if (x > xMin)
	    {
		if (x < yMin)
		{
		    iX++;
		    dist = ((d[iX] - distX) + (d[iY] - distY)) / 2;
		    d[iX] = dist;
		    iY += nleft - 2 - x;
		    r[x] += dist;
		    r[xMin] += dist;
		}
		else if (x > yMin)
		{
		    iX++;
		    iY++;
		    dist = ((d[iX] - distX) + (d[iY] - distY)) / 2;
		    d[iX] = dist;
		    r[x] += dist;
		    r[xMin] += dist;
		}
		else // if (x == yMin)
		{
		    iX++;
		}
	    }
	    else // if (x == xMin)
	    {
		iY += nleft - 2 - x;
	    }
	}

	/* Fill in the remaining gap (yMin row/column), by moving the first row
	 * into the gap.  The first row can be removed from the matrix in
	 * constant time, whereas collapsing the gap directly would require a
	 * series of memmove() calls, and leaving the gap would result in
	 * increased cache misses. */
	for (x = 1,
		 iX = CxpTreeNjXy2i(nleft, 0, x),
		 iY = CxpTreeNjXy2i(nleft, x, yMin);
	     x < nleft;
	     x++)
	{
	    if (x < yMin)
	    {
		d[iY] = d[iX];
		iY += nleft - 2 - x;
	    }
	    else if (x > yMin)
	    {
		iY++;
		d[iY] = d[iX];
	    }
	    iX++;
	}
	/* Fill in the gap in r, rScaled, and nodes. */
	r[yMin] = r[0];
	rScaled[yMin] = rScaled[0];
	nodes[yMin] = nodes[0];

	/* Move pointers forward, which removes the first row. */
	d = &d[nleft - 1];
	r = &r[1];
	rScaled = &rScaled[1];
	nodes = &nodes[1];
    }

    /* Join the remaining two nodes. */
#ifdef CxmTreeNjVerbose
    {
	PyObject *result;

	fprintf(stderr, "Join last two nodes:");

	result = CxNodeTaxonNumGet(nodes[0]);
	if (result != Py_None)
	{
	    fprintf(stderr, "%ld %p", PyInt_AsLong(result), nodes[0]);
	}
	else
	{
	    fprintf(stderr, "internal %p", nodes[0]);
	}

	result = CxNodeTaxonNumGet(nodes[1]);
	if (result != Py_None)
	{
	    fprintf(stderr, "and %ld %p\n", PyInt_AsLong(result), nodes[1]);
	}
	else
	{
	    fprintf(stderr, "and internal %p\n", nodes[1]);
	}
    }
#endif
    edge = CxEdgeNew(aTree);
    CxEdgeAttach(edge, nodes[0], nodes[1]);
    CxEdgeLengthSet(edge, d[0]);

    /* Set the tree base. */
    CxTreeBaseSet(aTree, nodes[0]);

    retval = false;
    RETURN:
    /* Clean up. */
    CxmFree(nodesOrig);
    CxmFree(rScaledOrig);
    CxmFree(rOrig);
    CxmFree(dOrig);
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
