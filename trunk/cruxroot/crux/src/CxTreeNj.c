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

typedef struct cw_tree_njd_s cw_tree_njd_t;
typedef struct cw_tree_njr_s cw_tree_njr_t;

/* Used by CxTree_nj() to represent a distance matrix cell. */
struct cw_tree_njd_s
{
    double dist;  /* Distance. */
    double trans; /* Transformed distance. */
};

/* Used by CxTree_nj(). */
struct cw_tree_njr_s
{
    double r;
    double r_scaled;   /* r/(m-2)). */
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
tr_p_nj_xy2i(uint32_t a_n, uint32_t a_x, uint32_t a_y)
{
    if (a_x > a_y)
    {
	uint32_t t;

	t = a_x;
	a_x = a_y;
	a_y = t;
    }

    return a_n * a_x + a_y - (((a_x + 3) * a_x) / 2) - 1;
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
//#define XXX_NJ_VERBOSE
static void
tree_p_nj(CxtCxtTreeObject *a_tree, double *a_distances, uint32_t a_ntaxa)
{
    cw_tree_njd_t *d, *d_prev, *t_d; /* Distance matrix. */
    cw_tree_njr_t *r; /* Distance sums. */
    uint32_t nleft, i, x, y, i_min, x_min, y_min, x_new, y_new;
    double dist_x, dist_y;
    CxtCxtNodeObject *node;
    CxtCxtEdgeObject *edge_x, *edge_y, *edge;

    cxmCheckPtr(a_distances);
    cxmAssert(a_ntaxa > 1);

    /* Allocate an array that is large enough to hold the distances. */
    d = (cw_tree_njd_t *) CxmMalloc(sizeof(cw_tree_njd_t)
				    * (tr_p_nj_xy2i(a_ntaxa,
						    a_ntaxa - 2,
						    a_ntaxa - 1)
				       + 1));
    /* d_prev can be smaller, since it won't be used until the matrix is shrunk
     * once. */
    d_prev = (cw_tree_njd_t *) CxmMalloc(sizeof(cw_tree_njd_t)
					 * (tr_p_nj_xy2i(a_ntaxa,
							 a_ntaxa - 3,
							 a_ntaxa - 2)
					    + 1));

    /* Initialize untransformed distances. */
    for (x = i = 0; x < a_ntaxa; x++)
    {
	for (y = x + 1; y < a_ntaxa; y++)
	{
	    d[i].dist = a_distances[x * a_ntaxa + y];

	    i++;
	}
    }

    /* Allocate an array that is large enough to hold all the distance sums. */
    r = (cw_tree_njr_t *) CxmMalloc(sizeof(cw_tree_njr_t) * a_ntaxa);

    /* Create a node for each taxon in the matrix. */
    for (i = 0; i < a_ntaxa; i++)
    {
	r[i].node = CxNodeNew(a_tree);
	CxNodeTaxonNumSetCargs(r[i].node, i);
    }

    /* Iteratitively join two nodes in the matrix, until only two are left. */
    for (nleft = a_ntaxa; nleft > 2; nleft--)
    {
#ifdef XXX_NJ_VERBOSE
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
	    r[i].r_scaled = r[i].r / (nleft - 2);
	}

	/* Calculate trans (transformed distance) for each pairwise distance.
	 * Keep track of the minimum trans, so that the corresponding nodes can
	 * be joined.  Ties are broken arbitrarily (the first minimum found is
	 * used). */
	for (x = i = 0, i_min = x_min = 0, y_min = 1; x < nleft; x++)
	{
	    for (y = x + 1; y < nleft; y++)
	    {
		d[i].trans = d[i].dist - ((r[x].r_scaled + r[y].r_scaled));

		if (d[i].trans < d[i_min].trans)
		{
		    i_min = i;
		    x_min = x;
		    y_min = y;
		}

		i++;
	    }
	}

#ifdef XXX_NJ_VERBOSE
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
	    fprintf(stderr, " || %8.4f\n", r[x].r_scaled);

	    fprintf(stderr, "\n");
	}
#endif
	

	/* Join the nodes with the minimum transformed distance. */
	node = CxNodeNew(a_tree);
#ifdef XXX_NJ_VERBOSE
	fprintf(stderr, "New node %ld %p\n",
		PyInt_AsLong(CxNodeTaxonNumGet(node)), node);
#endif
	edge_x = CxEdgeNew(a_tree);
	CxEdgeAttachCargs(edge_x, node, r[x_min].node);
	dist_x = (d[i_min].dist + r[x_min].r_scaled - r[y_min].r_scaled) / 2;
	CxEdgeLengthSetCargs(edge_x, dist_x);
#ifdef XXX_NJ_VERBOSE
	fprintf(stderr, "  Join node %ld %p (len %.4f)\n",
		PyInt_AsLong(CxNodeTaxonNumGet(r[x_min].node)), r[x_min].node,
		dist_x);
#endif

	edge_y = CxEdgeNew(a_tree);
	CxEdgeAttachCargs(edge_y, node, r[y_min].node);
	dist_y = d[i_min].dist - dist_x;
	CxEdgeLengthSetCargs(edge_y, dist_y);
#ifdef XXX_NJ_VERBOSE
	fprintf(stderr, "  Join node %ld %p (len %.4f)\n",
		PyInt_AsLong(CxNodeTaxonNumGet(r[y_min].node)), r[y_min].node,
		dist_y);
#endif

	/* Swap to new matrix. */
	t_d = d;
	d = d_prev;
	d_prev = t_d;

	/* Create compacted matrix. */
	for (x = i = 0, x_new = 0; x < nleft; x++)
	{
	    if (x != x_min && x != y_min)
	    {
		for (y = x + 1, y_new = x_new + 1; y < nleft; y++)
		{
		    if (y != x_min && y != y_min)
		    {
			d[tr_p_nj_xy2i(nleft - 1, x_new, y_new)].dist
			    = d_prev[i].dist;

			y_new++;
		    }

		    i++;
		}

		x_new++;
	    }
	    else
	    {
		i += (nleft - (x + 1));
	    }
	}

	/* Calculate distances to new node. */
	for (x = 0, x_new = 0; x < nleft; x++)
	{
	    if (x != x_min && x != y_min)
	    {
		d[tr_p_nj_xy2i(nleft - 1, x_new, nleft - 2)].dist
		    = ((d_prev[tr_p_nj_xy2i(nleft, x, x_min)].dist - dist_x)
		       + (d_prev[tr_p_nj_xy2i(nleft, x, y_min)].dist - dist_y)
		       ) / 2;

		x_new++;
	    }
	}

	/* Compact and update r. */
	// XXX Do this in a single pass?
	memmove(&r[y_min], &r[y_min + 1],
		sizeof(cw_tree_njr_t) * (nleft - y_min - 1));
	memmove(&r[x_min], &r[x_min + 1],
		sizeof(cw_tree_njr_t) * (nleft - x_min - 1));
	r[nleft - 2].node = node;
    }

    /* Join the remaining two nodes. */
#ifdef XXX_NJ_VERBOSE
    fprintf(stderr, "Join last two nodes: %ld %p and %ld %p\n",
	    PyInt_AsLong(CxNodeTaxonNumGet(r[0].node)), r[0].node,
	    PyInt_AsLong(CxNodeTaxonNumGet(r[1].node)), r[1].node);
#endif
    edge = CxEdgeNew(a_tree);
    CxEdgeAttachCargs(edge, r[0].node, r[1].node);
    CxEdgeLengthSetCargs(edge, d[0].dist);

    /* Set the tree base. */
    CxTreeBaseSetCargs(a_tree, r[0].node);

    /* Clean up. */
    CxmFree(d);
    CxmFree(d_prev);
    CxmFree(r);
}

PyObject *
CxTree_nj(CxtCxtTreeObject *self, PyObject *args)
{
    PyObject *retval, *dist_list, *tobj;
    double *distances;
    uint32_t ntaxa, nelms, i;
    bool okay;

    if (PyArg_ParseTuple(args, "O!", &PyList_Type, &dist_list) == 0)
    {
	retval = NULL;
	goto RETURN;
    }

    nelms = PyList_Size(dist_list);
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
	    tobj = PyList_GetItem(dist_list, i);
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
	    CxtTrNode old_tr_node, tr_node;
	    CxtCxtNodeObject *node;

	    old_tr_node = CxTrBaseGet(self->tr);

	    /* Neighbor-join. */
	    tree_p_nj(self, distances, ntaxa);

	    /* Reference new base. */
	    tr_node = CxTrBaseGet(self->tr);
	    if (tr_node != CxmTrNodeNone)
	    {
		node = (CxtCxtNodeObject *) CxTrNodeAuxGet(self->tr, tr_node);
		Py_INCREF(node);
	    }

	    /* Decref old base. */
	    if (old_tr_node != CxmTrNodeNone)
	    {
		node = (CxtCxtNodeObject *) CxTrNodeAuxGet(self->tr, old_tr_node);
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
