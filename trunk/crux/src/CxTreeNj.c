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
 * This file implements the neigbhor joining (NJ) and relaxed neighbor joining
 * (RNJ) algorithms.
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
 * Throughout the code, the matrix is stored as an upper-triangle symmetric
 * matrix, with additional bookkeeping, as necessary for RNJ and NJ.  For
 * example (n is current matrix size):
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
 * The key to RNJ's performance is the way in which clustering (node joining)
 * decisions are made.  Rather than calculating the transformed distances for
 * all possible node pairings and joining those two nodes as in NJ, RNJ
 * iteratively checks to see if various possible pairings of nodes are legal,
 * according to the constraints of the neighbor joining algorithm.  For example,
 * consider the joining of nodes 4 and 6:
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
 * As long as:
 *   1) the transformed distance for (4,6), denoted by **, is less than or equal
 *      to the transformed distances for the matrix elements marked by XX or YY,
 *      and
 *   2) joining rows 4 and 6 does not change the additivity of distances for the
 *      nodes still represented in the distance matrix (only applicable when
 *      dealing with additive distance matrices),
 * then joining nodes 4 and 6 poses no correctness problems for the neighbor
 * joining algorithm.  This implementation searches for such clusterings in an
 * efficient manner.
 *
 * It is important to note that in the worst case, this implementation has
 * O(n^3) performance.  However, worst case performance requires that the tree
 * be very long, with a particular pattern of branch lengths, and that the taxa
 * be inserted into the matrix in a particular order.  As such, worst case
 * performance almost never occurs, and typical runtime is approximately
 * proportional to (n^2 lg n).
 *
 ******************************************************************************/

#include "../include/_cruxmodule.h"

#include <float.h>

//#define CxmTreeNjDump

#ifdef CxmTreeNjDump
static void
CxpTreeNjDump(float *aD, float *aR, float *aRScaled, CxtNodeObject **aNodes,
	      long aNleft)
{
    float *dElm;
    long x, y;
    uint32_t taxonNum;
     
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
	taxonNum = CxNodeTaxonNumGet(aNodes[x]);
	if (taxonNum != CxmTrNodeTaxonNone)
	{
	    fprintf(stderr, " || n %u\n", taxonNum);
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
    Py_DECREF(aNodes[aXMin]);
    iMin = CxpTreeNjXy2i(aNleft, aXMin, aYMin);
    distX = (aD[iMin] + aRScaled[aXMin] - aRScaled[aYMin]) / 2;
    CxEdgeLengthSet(edgeX, distX);
    Py_DECREF(edgeX);

    edgeY = CxEdgeNew(aTree);
    CxEdgeAttach(edgeY, node, aNodes[aYMin]);
    Py_DECREF(aNodes[aYMin]);
    distY = aD[iMin] - distX;
    CxEdgeLengthSet(edgeY, distY);
    Py_DECREF(edgeY);

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
CxpTreeNjCompact(float *aD, float *aR, CxtNodeObject **aNodes, long aNleft,
		 long aXMin, long aYMin, CxtNodeObject *aNode,
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

    /* Fill in the gap in r, and nodes.  rScaled is re-calculated from scratch,
     * so there is no need to touch it here. */
    aR[aYMin] = aR[0];
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
    Py_DECREF(edge);
    Py_DECREF(aNodes[1]);

    return aNodes[0];
}

/* Compare two distances, and consider them equal if they are close enough. */
CxmpInline int
CxpTreeNjDistCompare(float aA, float aB)
{
    int retval;

    if (aB != 0.0)
    {
	if (fabs(1.0 - (aA / aB)) < FLT_EPSILON)
	{
	    retval = 0;
	}
	else if (aA < aB)
	{
	    retval = -1;
	}
	else
	{
	    retval = 1;
	}
    }
    else
    {
	if (aA <= -FLT_EPSILON)
	{
	    retval = -1;
	}
	else if (aA >= FLT_EPSILON)
	{
	    retval = 1;
	}
	else
	{
	    retval = 0;
	}
    }

    return retval;
}

CxmpInline void
CxpTreeNjRandomMinFind(float *aD, float *aRScaled, long aNleft, CxtMt *aMt,
		       long *rXMin, long *rYMin)
{
    float *dElm, transMin, transCur;
    long x, y, xMin, yMin, nmins;

    /* Calculate the transformed distance for each pairwise distance.  Keep
     * track of the minimum transformed distance, so that the corresponding
     * nodes can be joined.
     *
     * This is by far the most time-consuming portion of NJ.
     *
     * Use pointer arithmetic (dElm), rather than d[i].  This appears to reduce
     * register pressure on x86, and has a significant positive performance
     * impact. */
#ifdef CxmCcSilence
    xMin = yMin = 0;
#endif
    for (x = 0, nmins = 1, dElm = aD, transMin = HUGE_VAL; x < aNleft; x++)
    {
	for (y = x + 1; y < aNleft; y++)
	{
	    transCur = *dElm - (aRScaled[x] + aRScaled[y]);
	    dElm++;

	    /* Use CxpTreeNjDistCompare() in order to compare transformed
	     * distances, so that random selection is possible. */
	    switch (CxpTreeNjDistCompare(transCur, transMin))
	    {
		case -1:
		{
		    nmins = 1;
		    xMin = x;
		    yMin = y;
		    transMin = transCur;
		    break;
		}
		case 0:
		{
		    /* Choose such that all tied distances have an equal
		     * probability of being chosen. Otherwise, break ties
		     * arbitrarily (use the first minimum that was
		     * found). */
		    nmins++;
		    if (CxMtSint32RangeGet(aMt, nmins) == 0)
		    {
			xMin = x;
			yMin = y;
			transMin = transCur;
		    }
		    break;
		}
		case 1:
		{
		    break;
		}
		default:
		{
		    CxmNotReached();
		}
	    }
	}
    }

    *rXMin = xMin;
    *rYMin = yMin;
}

CxmpInline void
CxpTreeNjDeterministicMinFind(float *aD, float *aRScaled, long aNleft,
			      long *rXMin, long *rYMin)
{
    float *dElm, transMin, transCur;
    long x, y, xMin, yMin, nmins;

    /* Calculate the transformed distance for each pairwise distance.  Keep
     * track of the minimum transformed distance, so that the corresponding
     * nodes can be joined.
     *
     * This is by far the most time-consuming portion of NJ.
     *
     * Use pointer arithmetic (dElm), rather than d[i].  This appears to reduce
     * register pressure on x86, and has a significant positive performance
     * impact. */
#ifdef CxmCcSilence
    xMin = yMin = 0;
#endif
    for (x = 0, nmins = 1, dElm = aD, transMin = HUGE_VAL; x < aNleft; x++)
    {
	for (y = x + 1; y < aNleft; y++)
	{
	    transCur = *dElm - (aRScaled[x] + aRScaled[y]);
	    dElm++;

	    /* Since an arbitrary tie-breaking decision is being made
	     * anyway, don't bother using CxpTreeNjDistCompare() here. */
	    if (transCur < transMin)
	    {
		xMin = x;
		yMin = y;
		transMin = transCur;
	    }
	}
    }

    *rXMin = xMin;
    *rYMin = yMin;
}

CxmpInline long
CxpTreeNjRowAllMinFind(float *d, float *aRScaled, long aNleft,
		       long aX, CxtMt *aMt, float *rDist)
{
    long retval;
    float *dElm, dist, minDist;
    long y, nmins;

    minDist = HUGE_VAL;

    /* Find the minimum distance from the node on row aX to any other node that
     * comes before it in the matrix. */
    if (aX != 0)
    {
	for (y = 0,
		 dElm = &d[CxpTreeNjXy2i(aNleft, y, aX)];
	     y < aX;
	     y++)
	{
	    dist = *dElm - (aRScaled[y] + aRScaled[aX]);
	    dElm += (aNleft - 2 - y);

	    switch (CxpTreeNjDistCompare(dist, minDist))
	    {
		case -1:
		{
		    nmins = 1;
		    minDist = dist;
		    retval = y;
		    break;
		}
		case 0:
		{
		    /* Choose y such that all tied distances have an equal
		     * probability of being chosen. */
		    nmins++;
		    if (CxMtSint32RangeGet(aMt, nmins) == 0)
		    {
			retval = y;
		    }
		    break;
		}
		case 1:
		{
		    break;
		}
		default:
		{
		    CxmNotReached();
		}
	    }
	}
	CxmAssert(minDist != HUGE_VAL);
    }

    /* Find the minimum distance from the node on row aX to any other node that
     * comes after it in the matrix. */
    if (aX < aNleft - 1)
    {
	for (y = aX + 1,
		 dElm = &d[CxpTreeNjXy2i(aNleft, aX, y)];
	     y < aNleft;
	     y++)
	{
	    dist = *dElm - (aRScaled[aX] + aRScaled[y]);
	    dElm++;

	    switch (CxpTreeNjDistCompare(dist, minDist))
	    {
		case -1:
		{
		    nmins = 1;
		    minDist = dist;
		    retval = y;
		    break;
		}
		case 0:
		{
		    /* Choose y such that all tied distances have an equal
		     * probability of being chosen. */
		    nmins++;
		    if (CxMtSint32RangeGet(aMt, nmins) == 0)
		    {
			retval = y;
		    }
		    break;
		}
		case 1:
		{
		    break;
		}
		default:
		{
		    CxmNotReached();
		}
	    }
	}
    }
    CxmAssert(minDist != HUGE_VAL);

    *rDist = minDist;
    return retval;
}

CxmpInline bool
CxpTreeNjRowAllMinOk(float *d, float *aRScaled, long aNleft, long aX,
		     float aDist)
{
    bool retval;
    float *dElm, dist;
    long y;

    /* Make sure that aDist is <= any transformed distance in the row portion of
     * row aX. */
    if (aX + 1 < aNleft)
    {
	for (y = aX + 1,
		 dElm = &d[CxpTreeNjXy2i(aNleft, aX, y)];
	     y < aNleft;
	     y++)
	{
	    dist = *dElm - (aRScaled[aX] + aRScaled[y]);
	    dElm++;

	    if (CxpTreeNjDistCompare(dist, aDist) == -1)
	    {
		retval = false;
		goto RETURN;
	    }
	}
    }

    /* Make sure that aDist is <= any transformed distance in the column portion
     * of row aX. */
    if (aX != 0)
    {
	for (y = 0,
		 dElm = &d[CxpTreeNjXy2i(aNleft, y, aX)];
	     y < aX;
	     y++)
	{
	    dist = *dElm - (aRScaled[y] + aRScaled[aX]);
	    dElm += (aNleft - 2 - y);

	    if (CxpTreeNjDistCompare(dist, aDist) == -1)
	    {
		retval = false;
		goto RETURN;
	    }
	}
    }

    retval = true;
    RETURN:
    return retval;
}

CxmpInline long
CxpTreeNjRowMinFind(float *d, float *aRScaled, long aNleft, long x)
{
    long retval
#ifdef CxmCcSilence
	= -1
#endif
	;
    long y;
    float *dElm, dist, minDist;

    /* Find the minimum distance from the node on row x to any other node that
     * comes after it in the matrix. */
    for (y = x + 1,
	     dElm = &d[CxpTreeNjXy2i(aNleft, x, y)],
	     minDist = HUGE_VAL;
	 y < aNleft;
	 y++)
    {
	dist = *dElm - (aRScaled[x] + aRScaled[y]);
	dElm++;

	/* Don't bother using CxpTreeNjDistCompare() here, since this function
	 * is only used for deterministic RNJ. */
	if (dist < minDist)
	{
	    minDist = dist;
	    retval = y;
	}
    }
    CxmAssert(minDist != HUGE_VAL);
    CxmAssert(retval != -1);

    return retval;
}

/* Make sure that clustering aA and aB would not change the distances between
 * nodes.  This must be done in order to make sure that we get the true tree, in
 * the case that the distance matrix corresponds to precisely one tree
 * (distances are additive).  If the distances are non-additive though, there is
 * no need to do this check. */
CxmpInline bool
CxpTreeNjPairClusterAdditive(float *aD, float *aRScaled, long aNleft,
			     long aA, long aB)
{
    bool retval;
    long iAB, iA, iB, x;
    float distA, distB, dist;

    /* Calculate distances from {aA,aB} to the new node. */
    iAB = CxpTreeNjXy2i(aNleft, aA, aB);
    distA = (aD[iAB] + aRScaled[aA] - aRScaled[aB]) / 2;
    distB = aD[iAB] - distA;

    /* Calculate distances to the new node, and make sure that they are
     * consistent with the current distances. */

    /* Iterate over the row portion of distances for aA and aB. */
    if (aB + 1 < aNleft)
    {
	for (x = aB + 1,
		 iA = CxpTreeNjXy2i(aNleft, aA, aB + 1),
		 iB = CxpTreeNjXy2i(aNleft, aB, aB + 1);
	     x < aNleft;
	     x++)
	{
	    dist = ((aD[iA] - distA) + (aD[iB] - distB)) / 2;
	    iA++;
	    iB++;

	    if (CxpTreeNjDistCompare(dist + distA, aD[CxpTreeNjXy2i(aNleft,
								    aA, x)])
		!= 0)
	    {
		retval = false;
		goto RETURN;
	    }

	    if (CxpTreeNjDistCompare(dist + distB, aD[CxpTreeNjXy2i(aNleft,
								    aB, x)])
		!= 0)
	    {
		retval = false;
		goto RETURN;
	    }
	}
    }

    /* Iterate over the first column portion of distances for aA and aB. */
    for (x = 0,
	     iA = aA - 1,
	     iB = aB - 1;
	 x < aA;
	 x++)
    {
	dist = ((aD[iA] - distA) + (aD[iB] - distB)) / 2;
	iA += aNleft - 2 - x;
	iB += aNleft - 2 - x;

	if (CxpTreeNjDistCompare(dist + distA, aD[CxpTreeNjXy2i(aNleft, x, aA)])
	    != 0)
	{
	    retval = false;
	    goto RETURN;
	}

	if (CxpTreeNjDistCompare(dist + distB, aD[CxpTreeNjXy2i(aNleft, x, aB)])
	    != 0)
	{
	    retval = false;
	    goto RETURN;
	}
    }

    /* (x == aA) */
    iB += aNleft - 2 - x;
    x++;

    /* Iterate over the first row portion of distances for aA, and the second
     * column portion of distances for aB. */
    for (;
	 x < aB;
	 x++)
    {
	iA++;
	dist = ((aD[iA] - distA) + (aD[iB] - distB)) / 2;
	iB += aNleft - 2 - x;

	if (CxpTreeNjDistCompare(dist + distA, aD[CxpTreeNjXy2i(aNleft, aA, x)])
	    != 0)
	{
	    retval = false;
	    goto RETURN;
	}

	if (CxpTreeNjDistCompare(dist + distB, aD[CxpTreeNjXy2i(aNleft, x, aB)])
	    != 0)
	{
	    retval = false;
	    goto RETURN;
	}
    }

    retval = true;
    RETURN:
    return retval;
}

/* Finish checking whether it is okay to cluster rows aA and aB;
 * CxpTreeNjRowMinFind() or CxpTreeNjRowAllMinFind() has already done some of
 * the work by the time this function is called.
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

    /* Iterate over the row portion of distances for aB.  Distances for aA were
     * already checked before this function was called. */
    if (aB < aNleft - 1)
    {
	for (x = aB + 1,
		 iB = CxpTreeNjXy2i(aNleft, aB, aB + 1);
	     x < aNleft;
	     x++)
	{
	    dist = aD[iB] - (aRScaled[x] + aRScaled[aB]);
	    /* Don't bother using CxpTreeNjDistCompare() here, since this
	     * function is only used for deterministic RNJ. */
	    if (dist < distAB)
	    {
		retval = false;
		goto RETURN;
	    }
	    iB++;
	}
    }

    /* Iterate over the first column portion of distances for aA and aB. */
    for (x = 0,
	     iA = aA - 1,
	     iB = aB - 1;
	 x < aA;
	 x++)
    {
	dist = aD[iA] - (aRScaled[x] + aRScaled[aA]);
	/* Don't bother using CxpTreeNjDistCompare() here, since this function
	 * is only used for deterministic RNJ. */
	if (dist < distAB)
	{
	    retval = false;
	    goto RETURN;
	}

	dist = aD[iB] - (aRScaled[x] + aRScaled[aB]);
	/* Don't bother using CxpTreeNjDistCompare() here, since this function
	 * is only used for deterministic RNJ. */
	if (dist < distAB)
	{
	    retval = false;
	    goto RETURN;
	}

	iA += aNleft - 2 - x;
	iB += aNleft - 2 - x;
    }

    /* (x == aA) */
    iB += aNleft - 2 - x;
    x++;

    /* Iterate over the second column portion of distances for aB.  Distances
     * for aA were already checked before this function was called. */
    for (;
	 x < aB;
	 x++)
    {
	dist = aD[iB] - (aRScaled[x] + aRScaled[aB]);
	/* Don't bother using CxpTreeNjDistCompare() here, since this function
	 * is only used for deterministic RNJ. */
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

/* Get a seed for the C-based PRNG from Python's PRNG. */
static bool
CxpTreeNjSeedGet(long *rSeed)
{
    bool retval;
    PyObject *globals, *locals, *result, *seedCode;

    globals = PyEval_GetGlobals();
    if (globals == NULL)
    {
	retval = true;
	goto RETURN;
    }

    locals = Py_BuildValue("{sl}", "max", LONG_MAX);
    if (locals == NULL)
    {
	retval = true;
	goto RETURN;
    }
	
    seedCode = Py_CompileString("\
import random\n\
seed = random.randint(0, max)\n\
",
				"<string>",
				Py_file_input);
    if (seedCode == NULL)
    {
	Py_DECREF(locals);
	retval = true;
	goto RETURN;
    }

    result = PyEval_EvalCode((PyCodeObject *) seedCode, globals, locals);
    if (result == NULL)
    {
	Py_DECREF(seedCode);
	Py_DECREF(locals);
	retval = true;
	goto RETURN;
    }
    Py_DECREF(result);
    Py_DECREF(seedCode);

    result = PyDict_GetItemString(locals, "seed");
    if (result == NULL)
    {
	Py_DECREF(locals);
	retval = true;
	goto RETURN;
    }
    *rSeed = PyInt_AsLong(result);
    Py_DECREF(locals);

    retval = false;
    RETURN:
    return retval;
}

/* Try all clusterings of two rows in the matrix, in a random order.  Do this in
 * as cache-friendly a manner as possible (keeping in mind that the matrix is
 * stored in row-major form).  This means:
 *
 * 1) For each randomly chosen row (x), find the row which is the closest (y),
 *    according to transformed distances, by calling CxpTreeNjRowAllMinFind().
 *
 * 2) If the additivity constraint is enabled, check whether clustering x and y
 *    would violate additivity, by calling CxpTreeNjPairClusterAdditive().
 *
 * 2) Check whether it is okay to cluster x and y, by calling
 *    CxpTreeNjAllMinOk(). */
static bool
CxpTreeNjRandomCluster(float **arD, float *aR, float *aRScaled,
		       CxtNodeObject ***arNodes, long aNleft,
		       CxtTreeObject *aTree,
		       bool aAdditive)
{
    bool retval;
    long randomRow, closestRow, x, y;
    float distX, distY;
    float *d = *arD;
    CxtNodeObject *node;
    CxtNodeObject **nodes = *arNodes;
    CxtMt mt;
    CxtRi ri;
    long seed;
    float dist;
    bool clustered;

    if (CxpTreeNjSeedGet(&seed))
    {
	retval = true;
	goto RETURN;
    }
    CxMtNew(&mt);
    CxMtUint32Seed(&mt, seed);
    CxRiNew(&ri);
    CxRiInit(&ri, aNleft);

    clustered = true;
    while (true)
    {
	if (clustered == false)
	{
	    aAdditive = false;
	    CxRiInit(&ri, aNleft);
	}
	clustered = false;
	while (CxRiIndGet(&ri) < CxRiNintsGet(&ri))
	{
	    /* Randomly sample a row, without replacement. */
	    randomRow = CxRiRandomGet(&ri, &mt);

	    /* Find a row that is closest to x. */
	    closestRow = CxpTreeNjRowAllMinFind(d, aRScaled, aNleft, randomRow,
						&mt, &dist);

	    if (randomRow < closestRow)
	    {
		x = randomRow;
		y = closestRow;
	    }
	    else
	    {
		x = closestRow;
		y = randomRow;
	    }

	    /* Make sure that no row is closer to y than x is. */
	    if ((aAdditive == false
		 || CxpTreeNjPairClusterAdditive(d, aRScaled, aNleft, x, y))
		&& CxpTreeNjRowAllMinOk(d, aRScaled, aNleft, closestRow, dist))
	    {
		clustered = true;
#ifdef CxmTreeNjDump
		CxpTreeNjDump(d, aR, aRScaled, nodes, aNleft);
#endif
		CxpTreeNjNodesJoin(d, aRScaled, nodes, aTree, aNleft, x, y,
				   &node, &distX, &distY);
		CxpTreeNjRSubtract(d, aR, aNleft, x, y);
		CxpTreeNjCompact(d, aR, nodes, aNleft, x, y, node,
				 distX, distY);
		CxpTreeNjDiscard(&d, &aR, &aRScaled, &nodes, aNleft);
		aNleft--;
		CxpTreeNjRScaledUpdate(aRScaled, aR, aNleft);

		/* Shrinking the matrix may have reduced it to the point
		 * that the enclosing loop will no longer function
		 * correctly.  Check this condition here, in order to reduce
		 * branch overhead for the case where no join is done. */
		if (aNleft == 2)
		{
		    goto OUT;
		}

		CxRiInit(&ri, aNleft);
	    }
	}
    }
    OUT:

    CxRiDelete(&ri);
    CxMtDelete(&mt);

#ifdef CxmTreeNjDump
    CxpTreeNjDump(d, aR, aRScaled, nodes, aNleft);
#endif

    *arD = d;
    *arNodes = nodes;

    retval = false;
    RETURN:
    return retval;
}

/* Iteratively try all clusterings of two rows in the matrix.  Do this in a
 * cache-friendly manner (keeping in mind that the matrix is stored in row-major
 * form).  This means:
 *
 * 1) For each row (x), find the row after it which is the closest (y),
 *    according to transformed distances, by calling CxpTreeNjRowMinFind().
 *    This operation scans the row portion of the distances for x, which is a
 *    fast operation.
 *
 * 2) If the additivity constraint is enabled, check whether clustering x and y
 *    would violate additivity, by calling CxpTreeNjPairClusterAdditive().
 *
 * 2) Check whether it is okay to cluster x and y, by calling
 *    CxpTreeNjPairClusterOk().
 *
 * 3) If x and y can be clustered, do so, then immediately try to cluster with
 *    x again (as long as collapsing the matrix didn't move row x). */
static void
CxpTreeNjDeterministicCluster(float **arD, float *aR, float *aRScaled,
			      CxtNodeObject ***arNodes, long aNleft,
			      CxtTreeObject *aTree,
			      bool aAdditive)
{
    long x, y;
    float distX, distY;
    float *d = *arD;
    CxtNodeObject *node;
    CxtNodeObject **nodes = *arNodes;
    bool clustered;

    clustered = true;
    while (true)
    {
	if (clustered == false)
	{
	    aAdditive = false;
	}
	clustered = false;
	for (x = 0; x < aNleft - 1;) /* y indexes one past x. */
	{
	    y = CxpTreeNjRowMinFind(d, aRScaled, aNleft, x);

	    if ((aAdditive == false
		 || CxpTreeNjPairClusterAdditive(d, aRScaled, aNleft, x, y))
		&& CxpTreeNjPairClusterOk(d, aRScaled, aNleft, x, y))
	    {
		clustered = true;
#ifdef CxmTreeNjDump
		CxpTreeNjDump(d, aR, aRScaled, nodes, aNleft);
#endif
		CxpTreeNjNodesJoin(d, aRScaled, nodes, aTree, aNleft, x, y,
				   &node, &distX, &distY);
		CxpTreeNjRSubtract(d, aR, aNleft, x, y);
		CxpTreeNjCompact(d, aR, nodes, aNleft, x, y, node,
				 distX, distY);
		CxpTreeNjDiscard(&d, &aR, &aRScaled, &nodes, aNleft);
		aNleft--;
		CxpTreeNjRScaledUpdate(aRScaled, aR, aNleft);

		/* Shrinking the matrix may have reduced it to the point
		 * that the enclosing loop will no longer function
		 * correctly.  Check this condition here, in order to reduce
		 * branch overhead for the case where no join is done. */
		if (aNleft == 2)
		{
		    goto OUT;
		}

		/* The indexing of the matrix is shifted as a result of
		 * having removed the first row.  Set x such that joining
		 * with this row is immediately tried again.
		 *
		 * Note that if x is 0, then the row is now at (y - 1); in
		 * that case, stay on row 0. */
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
    }
    OUT:
#ifdef CxmTreeNjDump
    CxpTreeNjDump(d, aR, aRScaled, nodes, aNleft);
#endif

    *arD = d;
    *arNodes = nodes;
}

/* Create a tree from a pairwise distance matrix, using the NJ algorithm. */
static bool
CxpTreeNj(CxtTreeObject *aTree, float *aD, long aNtaxa, bool aRandom)
{
    bool retval;
    float *rOrig, *r; /* Distance sums. */
    float *rScaledOrig, *rScaled; /* Scaled distance sums: r/(nleft-2)). */
    CxtNodeObject **nodesOrig, **nodes; /* Nodes associated with each row. */
    CxtNodeObject *node;
    long nleft, xMin, yMin;
    float distX, distY;
    CxtMt mt;
    long seed;

    CxmCheckPtr(aD);
    CxmAssert(aNtaxa > 1);

    if (aRandom)
    {
        if (CxpTreeNjSeedGet(&seed))
	{
	    retval = true;
	    goto RETURN;
	}
	CxMtNew(&mt);
	CxMtUint32Seed(&mt, seed);
    }

    /* Initialize distance matrix, r, rScaled, and nodes. */
    rOrig = r = CxpTreeNjRInit(aD, aNtaxa);
    rScaledOrig = rScaled = CxpTreeNjRScaledInit(aNtaxa);
    CxpTreeNjRScaledUpdate(rScaled, r, aNtaxa);
    nodesOrig = nodes = CxpTreeNjNodesInit(aTree, aNtaxa);

    /* Iteratitively join two nodes in the matrix, until only two are left. */
    for (nleft = aNtaxa; nleft > 2; nleft--)
    {
	/* Standard neighbor joining. */
	CxpTreeNjRScaledUpdate(rScaled, r, nleft);
#ifdef CxmTreeNjDump
	CxpTreeNjDump(aD, r, rScaled, nodes, nleft);
#endif
	if (aRandom)
	{
	    CxpTreeNjRandomMinFind(aD, rScaled, nleft, &mt, &xMin, &yMin);
	}
	else
	{
	    CxpTreeNjDeterministicMinFind(aD, rScaled, nleft, &xMin, &yMin);
	}
	CxpTreeNjNodesJoin(aD, rScaled, nodes, aTree, nleft, xMin, yMin,
			   &node, &distX, &distY);
	CxpTreeNjRSubtract(aD, r, nleft, xMin, yMin);
	CxpTreeNjCompact(aD, r, nodes, nleft, xMin, yMin, node, distX, distY);
	CxpTreeNjDiscard(&aD, &r, &rScaled, &nodes, nleft);
    }

    /* Join last two nodes. */
    node = CxpTreeNjFinalJoin(aD, nodes, aTree);

    /* Set the tree base. */
    CxTreeBaseSet(aTree, node);
    Py_DECREF(node);

    retval = false;
    /* Clean up. */
    CxmFree(nodesOrig);
    CxmFree(rScaledOrig);
    CxmFree(rOrig);
    if (aRandom)
    {
        CxMtDelete(&mt);
    }
    RETURN:
    return retval;
}

/* Create a tree from a pairwise distance matrix, using the RNJ algorithm. */
static bool
CxpTreeRnj(CxtTreeObject *aTree, float *aD, long aNtaxa, bool aAdditive,
	   bool aRandom)
{
    bool retval;
    float *rOrig, *r; /* Distance sums. */
    float *rScaledOrig, *rScaled; /* Scaled distance sums: r/(nleft-2)). */
    CxtNodeObject **nodesOrig, **nodes; /* Nodes associated with each row. */
    CxtNodeObject *node;

    CxmCheckPtr(aD);
    CxmAssert(aNtaxa > 1);

    /* Initialize distance matrix, r, rScaled, and nodes. */
    rOrig = r = CxpTreeNjRInit(aD, aNtaxa);
    rScaledOrig = rScaled = CxpTreeNjRScaledInit(aNtaxa);
    CxpTreeNjRScaledUpdate(rScaled, r, aNtaxa);
    nodesOrig = nodes = CxpTreeNjNodesInit(aTree, aNtaxa);

    if (aRandom)
    {
	/* Cluster randomly. */
	if (CxpTreeNjRandomCluster(&aD, r, rScaled, &nodes, aNtaxa, aTree,
				   aAdditive))
	{
	    retval = true;
	    goto RETURN;
	}
    }
    else
    {
	/* Cluster deterministically. */
	CxpTreeNjDeterministicCluster(&aD, r, rScaled, &nodes, aNtaxa, aTree,
				      aAdditive);
    }

    /* Join last two nodes. */
    node = CxpTreeNjFinalJoin(aD, nodes, aTree);

    /* Set the tree base. */
    CxTreeBaseSet(aTree, node);
    Py_DECREF(node);

    retval = false;
    RETURN:
    /* Clean up. */
    CxmFree(nodesOrig);
    CxmFree(rScaledOrig);
    CxmFree(rOrig);
    return retval;
}

PyObject *
CxTreeNj(CxtTreeObject *self, PyObject *args)
{
    PyObject *retval, *distMatrix;
    float *d;
    long ntaxa;
    int random;

    if (PyArg_ParseTuple(args, "Oi", &distMatrix, &random) == 0)
    {
	retval = NULL;
	goto RETURN;
    }

    Py_INCREF(Py_None);
    retval = Py_None;
    CxmXepBegin();
    CxmXepTry
    {
	CxDistMatrixUpperHandoff((CxtDistMatrixObject *) distMatrix,
				 &d, &ntaxa);

	if (ntaxa > 1)
	{
	    /* NJ. */
	    if (CxpTreeNj(self, d, ntaxa, (bool) random))
	    {
		Py_DECREF(retval);
		retval = NULL;
	    }
	    CxmFree(d);
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

PyObject *
CxTreeRnj(CxtTreeObject *self, PyObject *args)
{
    PyObject *retval, *distMatrix;
    float *d;
    long ntaxa;
    int additive, random;

    if (PyArg_ParseTuple(args, "Oii", &distMatrix, &additive, &random) == 0)
    {
	retval = NULL;
	goto RETURN;
    }

    Py_INCREF(Py_None);
    retval = Py_None;
    CxmXepBegin();
    CxmXepTry
    {
	CxDistMatrixUpperHandoff((CxtDistMatrixObject *) distMatrix,
				 &d, &ntaxa);

	if (ntaxa > 1)
	{
	    /* RNJ. */
	    if (CxpTreeRnj(self, d, ntaxa, (bool) additive, (bool) random))
	    {
		Py_DECREF(retval);
		retval = NULL;
	    }
	    CxmFree(d);
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
