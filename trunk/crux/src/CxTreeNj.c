/******************************************************************************
 *
 * <Copyright = jasone>
 * <License>
 *
 ******************************************************************************
 *
 * Version: Crux <Version = crux>
 *
 ******************************************************************************
 *
 * This file implements the neighbor joining algorithm, but implements some
 * heuristic optimizations that dramatically improve performance, with no loss
 * of algorithmic correctness.
 *
 * Matrices are stored as upper-half matrices, which means that addressing is a
 * bit tricky.  The following matrix shows how addressing logically works, as
 * well as the order in which matrix elements are stored in memory:
 *   
 *   x: 0  1  2  3  4  5  6  7  8
 *     --+--+--+--+--+--+--+--+--+ y:
 *       | 0| 1| 2| 3| 4| 5| 6| 7| 0
 *       +--+--+--+--+--+--+--+--+
 *          | 8| 9|10|11|12|13|14| 1
 *          +--+--+--+--+--+--+--+
 *             |15|16|17|18|19|20| 2
 *             +--+--+--+--+--+--+
 *                |21|22|23|24|25| 3
 *                +--+--+--+--+--+
 *                   |26|27|28|29| 4
 *                   +--+--+--+--+
 *                      |30|31|32| 5
 *                      +--+--+--+
 *                         |33|34| 6
 *                         +--+--+
 *                            |35| 7
 *                            +--+
 *                               | 8
 *
 * The following formula can be used to convert from (x,y) coordinates to array
 * offsets:
 *
 *   n : Number of nodes currently in the matrix.
 *   x : Row.
 *   y : Column.
 *
 *                        2
 *                       x  + 3x
 *   f(n,x,y) = nx + y - ------- - 1
 *                          2
 *
 ******************************************************************************
 *
 * Since neighbor joining involves repeatedly joining two nodes and removing a
 * row from the matrix, the performance of matrix collapsing is important.
 * Additionally, it is important to keep the matrix compactly stored in memory,
 * for cache locality reasons.  This implementation removes row x by moving row
 * 0 into its place, then discarding row 0 of the array.  This has the effects
 * of 1) re-ordering rows, and 2) shifting row addresses, so care is necessary
 * in code that both iterates over rows and collapses the matrix.
 *
 ******************************************************************************
 *
 * The key to this implementation's performance is the way in which clustering
 * (node joining) decisions are made.  Rather than calculating the transformed
 * distances for all possible node pairings, and joining those two nodes, this
 * implementation iteratively checks to see if various possible pairings of
 * nodes are legal, according to the constraints of the neighbor joining
 * algorithm.  For example, consider the joining of nodes 4 and 6:
 *
 *   x: 0  1  2  3  4  5  6  7  8
 *     --+--+--+--+--+--+--+--+--+ y:
 *       |  |  |  |XX|  |YY|  |  | 0
 *       +--+--+--+--+--+--+--+--+
 *          |  |  |XX|  |YY|  |  | 1
 *          +--+--+--+--+--+--+--+
 *             |  |XX|  |YY|  |  | 2
 *             +--+--+--+--+--+--+
 *                |XX|  |YY|  |  | 3
 *                +--+--+--+--+--+
 *                   |XX|**|XX|XX| 4
 *                   +--+--+--+--+
 *                      |YY|  |  | 5
 *                      +--+--+--+
 *                         |YY|YY| 6
 *                         +--+--+
 *                            |  | 7
 *                            +--+
 *                               | 8
 *
 * As long as the transformed distance for (4,6), denoted by **, is less than or
 * equal to the transformed distances for the matrix elements marked by XX or
 * YY, then joining nodes 4 and 6 poses no correctness problems for the neighbor
 * joining algorithm.  This implementation searches for such clusterings in an
 * efficient manner.
 *
 * It is important to note that in the worst case, this implementation has
 * O(n^3) performance.  However, worst case performance requires that the tree
 * be very long, with a particular pattern of branch lengths, and that the taxa
 * be inserted into the matrix in a particular order.  As such, worst case
 * performance almost never occurs.
 *
 ******************************************************************************/

#include "../include/_cruxmodule.h"

//#define CxmTreeNjRandomize
//#define CxmTreeNjVerbose
//#define CxmTreeNjDump

#ifdef CxmTreeNjDump
static void
CxpTreeNjDump(float *aD, float *aR, float *aRScaled, CxtNodeObject **aNodes,
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
	fprintf(stderr, " || %8.4f\n", aR[x]);

	fprintf(stderr, "%*s", (int) x * 9 + (!!x), " ");
	for (y -= aNleft - (x + 1), dElm -= aNleft - (x + 1);
	     y < aNleft;
	     y++)
	{
	    fprintf(stderr, " %8.4f", *dElm - (aRScaled[x] + aRScaled[y]));
	    dElm++;
	}
	fprintf(stderr, " || %8.4f\n", aRScaled[x]);

	fprintf(stderr, "%*s", (int) x * 9 + (!!x), " ");
	for (y -= aNleft - (x + 1);
	     y < aNleft;
	     y++)
	{
	    fprintf(stderr, " %8s", "");
	}
	result = CxNodeTaxonNumGet(aNodes[x]);
	if (result != Py_None)
	{
	    fprintf(stderr, " || n %ld\n",
		    PyInt_AsLong(result), aNodes[x]);
	}
	else
	{
	    fprintf(stderr, " || i %p\n", aNodes[x]);
	}
    }
}
#endif

/* Convert from row/column matrix coordinates to array offsets. */
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

    // XXX Use dElm in this function and others that do not.

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

/* Finish checking whether it is okay to cluster rows aA and aB;
 * CxpTreeNjCluster() has already done some of the work by the time this
 * function is called.
 *
 * Two nodes, aA and aB, can be clustered if the transformed distance between
 * them is less than or equal to the transformed distances from aA or aB to any
 * other node.
 */
CxmpInline bool
CxpTreeNjPairClusterOk(float *aD, float *aRScaled, long aNleft,
		       long aA, long aB)
{
    bool retval;
    long x, iA, iB;
    float distAB, dist;

    CxmAssert(aA < aB);

    /* Calculate the transformed distance between aA and aB. */
    distAB = aD[CxpTreeNjXy2i(aNleft, aA, aB)] - (aRScaled[aA] + aRScaled[aB]);

    /* Iterate over the row-major portion of distances for aB.  Distances for aA
     * were already checked in CxpTreeNjCluster(). */
    if (aB < aNleft - 1)
    {
	for (x = aB + 1,
		 iB = CxpTreeNjXy2i(aNleft, aB, aB + 1);
	     x < aNleft;
	     x++)
	{
	    dist = aD[iB] - (aRScaled[x] + aRScaled[aB]);
	    if (dist < distAB)
	    {
		retval = false;
		goto RETURN;
	    }
	    iB++;
	}
    }

    /* Iterate over the first column-major portion of distances for aA and
     * aB. */
    for (x = 0,
	     iA = aA - 1,
	     iB = aB - 1;
	 x < aA;
	 x++)
    {
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

    /* Iterate over the second column-major portion of distances for aB.
     * Distances for aA were already checked in CxpTreeNjCluster(). */
    for (;
	 x < aB;
	 x++)
    {
	dist = aD[iB] - (aRScaled[x] + aRScaled[aB]);
	if (dist < distAB)
	{
	    retval = false;
	    goto RETURN;
	}
	iB += aNleft - 2 - x;
    }

    retval = true;
    RETURN:
    return retval;
}

/* Compare two distances, and consider them equal if they are close enough. */
CxmpInline bool
CxpTreeNjDistEq(float aA, float aB)
{
    bool retval;
    float ratio;
    // XXX What is a reasonable value for this?  Rounding error appears to start
    // showing up at 1.0e-7.
#define CxmTreeNjMaxDiff 1.0e-6

    ratio = aA / aB;
    if (ratio < (1.0 - CxmTreeNjMaxDiff) || ratio > (1.0 + CxmTreeNjMaxDiff))
    {
//	fprintf(stderr, "%.20e != %.20e (ratio: %.10e)\n", aA, aB, ratio);
	retval = false;
	goto RETURN;
    }

    retval = true;
    RETURN:
    return retval;
}

// XXX It may be enough to check a single distance.  If so, vastly simplify
// this function.

/* Make sure that clustering aA and aB would not change the distances between
 * nodes.  This must be done in order to make sure that we get the true tree, in
 * the case that the distance matrix corresponds to precisely one tree.
 * If the distance matrix is inconsistent though, there is no need to do this
 * check. */
CxmpInline bool
CxpTreeNjPairClusterExact(float *aD, float *aRScaled, long aNleft,
			  long aA, long aB)
{
    bool retval;
    long iAB, x, iA, iB;
    float distA, distB, dist;

    /* Calculate distances from {aA,aB} to the new node. */
    iAB = CxpTreeNjXy2i(aNleft, aA, aB);
    distA = (aD[iAB] + aRScaled[aA] - aRScaled[aB]) / 2;
    distB = aD[iAB] - distA;

    // XXX Reverse order of following three loops, for improved cache
    // performance.

    /* Calculate distances to the new node, and make sure that they are
     * consistent with the current distances. */
    for (x = 0,
	     iA = aA - 1,
	     iB = aB - 1;
	 x < aA;
	 x++)
    {
	dist = ((aD[iA] - distA) + (aD[iB] - distB)) / 2;
	iA += aNleft - 2 - x;
	iB += aNleft - 2 - x;

	if (CxpTreeNjDistEq(dist + distA, aD[CxpTreeNjXy2i(aNleft, x, aA)])
	    == false)
	{
// 	    fprintf(stderr, "(%ld,%ld): Unequal distances (%ld, %ld)\n",
// 		    aA, aB, x, aA);
	    retval = false;
	    goto RETURN;
	}

	if (CxpTreeNjDistEq(dist + distB, aD[CxpTreeNjXy2i(aNleft, x, aB)])
	    == false)
	{
// 	    fprintf(stderr, "(%ld,%ld): Unequal distances (%ld, %ld)\n",
// 		    aA, aB, x, aA);
	    retval = false;
	    goto RETURN;
	}
    }

    /* (x == aA) */
    iB += aNleft - 2 - x;
    x++;

    for (;
	 x < aB;
	 x++)
    {
	iA++;
	dist = ((aD[iA] - distA) + (aD[iB] - distB)) / 2;
	iB += aNleft - 2 - x;

	if (CxpTreeNjDistEq(dist + distA, aD[CxpTreeNjXy2i(aNleft, aA, x)])
	    == false)
	{
// 	    fprintf(stderr, "(%ld,%ld): Unequal distances (%ld, %ld)\n",
// 		    aA, aB, x, aA);
	    retval = false;
	    goto RETURN;
	}

	if (CxpTreeNjDistEq(dist + distB, aD[CxpTreeNjXy2i(aNleft, x, aB)])
	    == false)
	{
// 	    fprintf(stderr, "(%ld,%ld): Unequal distances (%ld, %ld)\n",
// 		    aA, aB, x, aA);
	    retval = false;
	    goto RETURN;
	}
    }

    /* (x == aB) */
    iA++;
    x++;

    for (;
	 x < aNleft;
	 x++)
    {
	iA++;
	iB++;
	dist = ((aD[iA] - distA) + (aD[iB] - distB)) / 2;

	if (CxpTreeNjDistEq(dist + distA, aD[CxpTreeNjXy2i(aNleft, aA, x)])
	    == false)
	{
// 	    fprintf(stderr, "(%ld,%ld): Unequal distances (%ld, %ld)\n",
// 		    aA, aB, x, aA);
	    retval = false;
	    goto RETURN;
	}

	if (CxpTreeNjDistEq(dist + distB, aD[CxpTreeNjXy2i(aNleft, aB, x)])
	    == false)
	{
// 	    fprintf(stderr, "(%ld,%ld): Unequal distances (%ld, %ld)\n",
// 		    aA, aB, x, aA);
	    retval = false;
	    goto RETURN;
	}
    }

    retval = true;
    RETURN:
    return retval;
}

/* Iteratively try all clusterings of two rows in the matrix.  Do this in a
 * cache-friendly manner (keeping in mind that the matrix is stored in row-major
 * form).  This means:
 *
 * 1) For each row (x), find the row after it which is the closest (min),
 *    according to transformed distances.  This operation scans the row-major
 *    portion of the distances for x, which is a fast operation.
 *
 * 2) Check whether it is okay to cluster x and min, by calling
 *    CxpTreeNjPairClusterOk().
 *
 * 3) If x and min can be clustered, do so, then immediately try to cluster with
 *    x again (as long as collapsing the matrix didn't move row x). */
static bool
CxpTreeNjCluster(float **arD, float **arR, float **arRScaled,
		 CxtNodeObject ***arNodes, long *arNleft, CxtTreeObject *aTree,
		 bool aExact)
{
    bool retval = false;
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
	/* Find the minimum distance from the node on row x to any other node
	 * that comes after it in the matrix.  This has the effect of trying
	 * each node pairing only once. */
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

	if (CxpTreeNjPairClusterOk(d, rScaled, nleft, x, min)
	    && (aExact == false ||
		(CxpTreeNjPairClusterExact(d, rScaled, nleft, x, min))))
	{
	    retval = true;
	    // XXX Move randomization to matrix initialization, and expose it as
	    // an option.
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

		if (rand() & 1)
		{
		    x++;
		    continue;
		}
	    }
#endif
#ifdef CxmTreeNjDump
	    CxpTreeNjDump(d, r, rScaled, nodes, nleft);
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
             * removed the first row.  Set x such that joining with this row is
             * immediately tried again.  This isn't ideal, in that this only
             * tries to join the new node with nodes that come after it in the
             * matrix.  However, it probably isn't worth the cache miss penalty
             * of finding the node that is closest to this one (requires
             * iterating over a column).
	     *
	     * Note that if x is 0, then the row is now at (min - 1); in that
	     * case, stay on row 0. */
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
    if (aExact == false)
    {
	retval = false;
    }
    // XXX Remove.
//     if (retval == false)
//     {
// 	fprintf(stderr, "Next round will be in inexact mode\n");
//     }
//     else
//     {
// 	fprintf(stderr, "Next round will be in exact mode\n");
//     }
    return retval;
}

/* Create a tree from a pairwise distance matrix, using the neighbor-joining
 * algorithm.
 *
 * The matrix is stored as an upper-triangle symmetric matrix, with additional
 * bookkeeping, as necessary for the neighbor-joining algorithm.  For example (n
 * is current matrix size):
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
    bool exact;
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
    exact = true;
    CxpTreeNjRScaledUpdate(rScaled, r, nleft);
    while (nleft > 2)
    {
	exact = CxpTreeNjCluster(&d, &r, &rScaled, &nodes, &nleft, aTree,
				 exact);
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
