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

typedef struct CxsTreeNjd CxtTreeNjd;
typedef struct CxsTreeNjr CxtTreeNjr;

/* Used by CxTreeNj() to represent a distance matrix cell. */
struct CxsTreeNjd
{
    double dist;  /* Distance. */
    double trans; /* Transformed distance. */
};

/* Used by CxTreeNj(). */
struct CxsTreeNjr
{
    double r;
    double rScaled;   /* r/(m-2)). */
    CxtCxtNodeObject *node; /* Associated node. */
};

/* The following function can be used to convert from row/column matrix
 * coordinates to array offsets (n is ntaxa, and x < y) for neighbor-joining:
 *
 *                        2
 *                       x  + 3x
 *   f(n,x,y) = nx + y - ------- - 1
 *                          2
 */
CxmpInline uint32_t
CxpTreeNjXy2i(uint32_t aN, uint32_t aX, uint32_t aY)
{
    if (aX > aY)
    {
	uint32_t t;

	t = aX;
	aX = aY;
	aY = t;
    }

    return aN * aX + aY - (((aX + 3) * aX) / 2) - 1;
}

/*
 * Create a tree from a pairwise distance matrix, using the neighbor-joining
 * algorithm.
 *
 * The matrix is actually stored as an upper-triangle symmetric matrix, with
 * additional bookkeeping, as necessary for the neighbor-joining algorithm.
 * For example (m is current matrix size):
 *
 *                                                    |   r   ||
 *                                                    |  ---  ||
 *   |   A   |   B   |   C   |   D   |   E   ||   r   |  m-2  ||   
 *   +=======+=======+=======+=======+=======++=======+=======++===
 *           |   0   |   1   |   2   |   3   ||   6   |   2.0 || A 
 *           | -10.5 | -12.0 | -12.5 | -13.0 ||       |       ||   
 *           +-------+-------+-------+-------++-------+-------++---
 *                   |   4   |   5   |   6   ||  15   |   5.0 || B 
 *                   | -13.5 | -14.0 | -14.5 ||       |       ||   
 *                   +-------+-------+-------++-------+-------++---
 *                           |   7   |   8   ||  20   |   6.7 || C 
 *                           | -14.5 | -15.0 ||       |       ||   
 *                           +-------+-------++-------+-------++---
 *                                   |   9   ||  23   |   7.7 || D 
 *                                   | -15.5 ||       |       ||   
 *                                   +-------++-------+-------++---
 *                                           ||  26   |   8.7 || E 
 *                                           ||       |       ||   
 *                                           ++-------+-------++---
 *
 * is stored as:
 *
 *   d:
 *   /-------+-------+-------+-------+-------+-------+-------+-------+-------...
 *   |   0   |   1   |   2   |   3   |   4   |   5   |   6   |   7   |   8   
 *   | -10.5 | -12.0 | -12.5 | -13.0 | -13.5 | -14.0 | -14.5 | -14.5 | -15.0 
 *   \-------+-------+-------+-------+-------+-------+-------+-------+-------...
 *
 *   +-------\
 *   |   9   |
 *   | -15.5 |
 *   +-------/
 *
 *   r:
 *   /-------+-------+-------+-------+-------\
 *   |   6   |  15   |  20   |  23   |  26   |	
 *   |   2.0 |   5.0 |   6.7 |   7.7 |   8.7 |	
 *   \-------+-------+-------+-------+-------/
 */
//#define CxmTreeNjVerbose
static void
CxpTreeNj(CxtCxtTreeObject *aTree, double *aDistances, uint32_t aNtaxa)
{
    CxtTreeNjd *d, *dPrev, *tD; /* Distance matrix. */
    CxtTreeNjr *r; /* Distance sums. */
    uint32_t nleft, i, x, y, iMin, xMin, yMin, xNew, yNew;
    double distX, distY;
    CxtCxtNodeObject *node;
    CxtCxtEdgeObject *edgeX, *edgeY, *edge;

    cxmCheckPtr(aDistances);
    cxmAssert(aNtaxa > 1);

    /* Allocate an array that is large enough to hold the distances. */
    d = (CxtTreeNjd *) CxmMalloc(sizeof(CxtTreeNjd)
				    * (CxpTreeNjXy2i(aNtaxa,
						    aNtaxa - 2,
						    aNtaxa - 1)
				       + 1));
    /* dPrev can be smaller, since it won't be used until the matrix is shrunk
     * once. */
    dPrev = (CxtTreeNjd *) CxmMalloc(sizeof(CxtTreeNjd)
					 * (CxpTreeNjXy2i(aNtaxa,
							 aNtaxa - 3,
							 aNtaxa - 2)
					    + 1));

    /* Initialize untransformed distances. */
    for (x = i = 0; x < aNtaxa; x++)
    {
	for (y = x + 1; y < aNtaxa; y++)
	{
	    d[i].dist = aDistances[x * aNtaxa + y];

	    i++;
	}
    }

    /* Allocate an array that is large enough to hold all the distance sums. */
    r = (CxtTreeNjr *) CxmMalloc(sizeof(CxtTreeNjr) * aNtaxa);

    /* Create a node for each taxon in the matrix. */
    for (i = 0; i < aNtaxa; i++)
    {
	r[i].node = CxNodeNew(aTree);
	CxNodeTaxonNumSetCargs(r[i].node, i);
    }

    /* Iteratitively join two nodes in the matrix, until only two are left. */
    for (nleft = aNtaxa; nleft > 2; nleft--)
    {
#ifdef CxmTreeNjVerbose
	fprintf(stderr, "Size: %u\n", nleft);
#endif

	/* Calculate r (sum of distances to other nodes) and r/(nleft-2)
	 * for each node. */
	for (i = 0; i < nleft; i++)
	{
	    r[i].r = 0.0;
	}

	for (x = i = 0; x < nleft; x++)
	{
	    for (y = x + 1; y < nleft; y++)
	    {
		r[x].r += d[i].dist;
		r[y].r += d[i].dist;

		i++;
	    }
	}

	for (i = 0; i < nleft; i++)
	{
	    r[i].rScaled = r[i].r / (nleft - 2);
	}

	/* Calculate trans (transformed distance) for each pairwise distance.
	 * Keep track of the minimum trans, so that the corresponding nodes can
	 * be joined.  Ties are broken arbitrarily (the first minimum found is
	 * used). */
	for (x = i = 0, iMin = xMin = 0, yMin = 1; x < nleft; x++)
	{
	    for (y = x + 1; y < nleft; y++)
	    {
		d[i].trans = d[i].dist - ((r[x].rScaled + r[y].rScaled));

		if (d[i].trans < d[iMin].trans)
		{
		    iMin = i;
		    xMin = x;
		    yMin = y;
		}

		i++;
	    }
	}

#ifdef CxmTreeNjVerbose
	fprintf(stderr,
		"----------------------------------------"
		"----------------------------------------\n");
	for (x = i = 0; x < nleft; x++)
	{
	    fprintf(stderr, "%*s", x * 9 + (!!x), " ");
	    for (y = x + 1; y < nleft; y++)
	    {
		fprintf(stderr, " %8.4f", d[i].dist);

		i++;
	    }
	    fprintf(stderr, " || %8.4f (node %ld %p)\n", r[x].r,
		    PyInt_AsLong(CxNodeTaxonNumGet(r[x].node)), r[x].node);

	    fprintf(stderr, "%*s", x * 9 + (!!x), " ");
	    for (i -= (nleft - (x + 1)), y = x + 1; y < nleft; y++)
	    {
		fprintf(stderr, " %8.4f", d[i].trans);

		i++;
	    }
	    fprintf(stderr, " || %8.4f\n", r[x].rScaled);

	    fprintf(stderr, "\n");
	}
#endif
	

	/* Join the nodes with the minimum transformed distance. */
	node = CxNodeNew(aTree);
#ifdef CxmTreeNjVerbose
	fprintf(stderr, "New node %ld %p\n",
		PyInt_AsLong(CxNodeTaxonNumGet(node)), node);
#endif
	edgeX = CxEdgeNew(aTree);
	CxEdgeAttachCargs(edgeX, node, r[xMin].node);
	distX = (d[iMin].dist + r[xMin].rScaled - r[yMin].rScaled) / 2;
	CxEdgeLengthSetCargs(edgeX, distX);
#ifdef CxmTreeNjVerbose
	fprintf(stderr, "  Join node %ld %p (len %.4f)\n",
		PyInt_AsLong(CxNodeTaxonNumGet(r[xMin].node)), r[xMin].node,
		distX);
#endif

	edgeY = CxEdgeNew(aTree);
	CxEdgeAttachCargs(edgeY, node, r[yMin].node);
	distY = d[iMin].dist - distX;
	CxEdgeLengthSetCargs(edgeY, distY);
#ifdef CxmTreeNjVerbose
	fprintf(stderr, "  Join node %ld %p (len %.4f)\n",
		PyInt_AsLong(CxNodeTaxonNumGet(r[yMin].node)), r[yMin].node,
		distY);
#endif

	/* Swap to new matrix. */
	tD = d;
	d = dPrev;
	dPrev = tD;

	/* Create compacted matrix. */
	for (x = i = 0, xNew = 0; x < nleft; x++)
	{
	    if (x != xMin && x != yMin)
	    {
		for (y = x + 1, yNew = xNew + 1; y < nleft; y++)
		{
		    if (y != xMin && y != yMin)
		    {
			d[CxpTreeNjXy2i(nleft - 1, xNew, yNew)].dist
			    = dPrev[i].dist;

			yNew++;
		    }

		    i++;
		}

		xNew++;
	    }
	    else
	    {
		i += (nleft - (x + 1));
	    }
	}

	/* Calculate distances to new node. */
	for (x = 0, xNew = 0; x < nleft; x++)
	{
	    if (x != xMin && x != yMin)
	    {
		d[CxpTreeNjXy2i(nleft - 1, xNew, nleft - 2)].dist
		    = ((dPrev[CxpTreeNjXy2i(nleft, x, xMin)].dist - distX)
		       + (dPrev[CxpTreeNjXy2i(nleft, x, yMin)].dist - distY)
		       ) / 2;

		xNew++;
	    }
	}

	/* Compact and update r. */
	// XXX Do this in a single pass?
	memmove(&r[yMin], &r[yMin + 1],
		sizeof(CxtTreeNjr) * (nleft - yMin - 1));
	memmove(&r[xMin], &r[xMin + 1],
		sizeof(CxtTreeNjr) * (nleft - xMin - 1));
	r[nleft - 2].node = node;
    }

    /* Join the remaining two nodes. */
#ifdef CxmTreeNjVerbose
    fprintf(stderr, "Join last two nodes: %ld %p and %ld %p\n",
	    PyInt_AsLong(CxNodeTaxonNumGet(r[0].node)), r[0].node,
	    PyInt_AsLong(CxNodeTaxonNumGet(r[1].node)), r[1].node);
#endif
    edge = CxEdgeNew(aTree);
    CxEdgeAttachCargs(edge, r[0].node, r[1].node);
    CxEdgeLengthSetCargs(edge, d[0].dist);

    /* Set the tree base. */
    CxTreeBaseSetCargs(aTree, r[0].node);

    /* Clean up. */
    CxmFree(d);
    CxmFree(dPrev);
    CxmFree(r);
}

PyObject *
CxTreeNj(CxtCxtTreeObject *self, PyObject *args)
{
    PyObject *retval, *distList, *tobj;
    double *distances;
    uint32_t ntaxa, nelms, i;
    bool okay;

    if (PyArg_ParseTuple(args, "O!", &PyList_Type, &distList) == 0)
    {
	retval = NULL;
	goto RETURN;
    }

    nelms = PyList_Size(distList);
    ntaxa = sqrt(nelms);
    if (ntaxa * ntaxa != nelms)
    {
	Py_INCREF(PyExc_ValueError);
	retval = PyExc_ValueError;
	goto RETURN;
    }

    CxmXepBegin();
    xep_try
    {
	distances = (double *) CxmMalloc(sizeof(double) * nelms);

	okay = true;
	for (i = 0; i < nelms; i++)
	{
	    tobj = PyList_GetItem(distList, i);
	    if (PyFloat_Check(tobj))
	    {
		distances[i] = PyFloat_AsDouble(tobj);
	    }
	    else if (PyInt_Check(tobj))
	    {
		distances[i] = PyInt_AsLong(tobj);
	    }
	    else
	    {
		Py_INCREF(PyExc_ValueError);
		retval = PyExc_ValueError;
		okay = false;
	    }
	}

	if (okay)
	{
	    CxtTrNode oldTrNode, trNode;
	    CxtCxtNodeObject *node;

	    oldTrNode = CxTrBaseGet(self->tr);

	    /* Neighbor-join. */
	    CxpTreeNj(self, distances, ntaxa);

	    /* Reference new base. */
	    trNode = CxTrBaseGet(self->tr);
	    if (trNode != CxmTrNodeNone)
	    {
		node = (CxtCxtNodeObject *) CxTrNodeAuxGet(self->tr, trNode);
		Py_INCREF(node);
	    }

	    /* Decref old base. */
	    if (oldTrNode != CxmTrNodeNone)
	    {
		node = (CxtCxtNodeObject *) CxTrNodeAuxGet(self->tr, oldTrNode);
		Py_DECREF(node);
	    }
	}

	CxmFree(distances);
    }
    CxmXepCatch(CxmXepOOM)
    {
	CxmXepHandled();
	retval = PyErr_NoMemory();
    }
    CxmXepEnd();

    Py_INCREF(Py_None);
    retval = Py_None;
    RETURN:
    return retval;
}
