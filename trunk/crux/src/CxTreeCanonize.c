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

// Used for canonizing trees.
struct CxsTreeCanonize
{
    CxtRingObject *ring;
    unsigned minTaxon;
};

// Comparison function that is passed to qsort().
static int
CxpTreeCanonizeCompare(const void *aA, const void *aB)
{
    const struct CxsTreeCanonize *a = (const struct CxsTreeCanonize *) aA;
    const struct CxsTreeCanonize *b = (const struct CxsTreeCanonize *) aB;

    if (a->minTaxon < b->minTaxon)
    {
	return -1;
    }
    else
    {
	CxmAssert(a->minTaxon > b->minTaxon);
	return 1;
    }
}

// Convert a tree node to canonical form by re-ordering the ring such that
// subtrees are in increasing order of minimum taxon number contained.
static unsigned
CxpTreeCanonize(CxtTreeObject *aTree, CxtRingObject *aRing)
{
    unsigned rVal, degree;
    CxtNodeObject *node;

    // Get taxon number (an internal node has CxmTrNodeTaxonNone).
    node = CxRingNode(aRing);
    rVal = CxNodeTaxonNumGet(node);

    // Get the degree of the node that this ring is a part of.
    degree = CxNodeDegree(node);

    if (degree > 1)
    {
	unsigned i, minTaxon;
	CxtNodeObject *nodeOther;
	CxtEdgeObject *edge;
	CxtRingObject *ring;
	struct CxsTreeCanonize canonize[degree - 1];

	// Iteratively canonize subtrees, keeping track of the minimum taxon
	// number seen overall, as well as for each subtree.
	i = 0;
	rVal = CxmTrNodeTaxonNone;
	for (ring = CxRingNext(aRing); ring != aRing; ring = CxRingNext(ring))
	{
	    minTaxon = CxpTreeCanonize(aTree, CxRingOther(ring));
	    if (minTaxon < rVal)
	    {
		rVal = minTaxon;
	    }

	    canonize[i].ring = ring;
	    canonize[i].minTaxon = minTaxon;

	    i++;
	}
	CxmAssert(i == degree - 1);

	// Sort the subtrees.
	qsort(canonize, degree - 1, sizeof(struct CxsTreeCanonize),
	      CxpTreeCanonizeCompare);

	// Set the beginning of the ring to aRing.  This makes it easier for
	// external code to traverse a tree in canonical order.
	CxNodeRingSet(node, aRing);

	// Detach and re-attach all edges, in order.  This code assumes that
	// attaching inserts the edge at the end of the ring.  The first element
	// can be skipped, since the removal/re-insertion of all other elements
	// eventually leaves the first element in the proper location.
	Py_INCREF(node);
	for (i = 1; i < (degree - 1); i++)
	{
	    nodeOther = CxRingNode(CxRingOther(canonize[i].ring));
	    edge = CxRingEdge(canonize[i].ring);

	    Py_INCREF(edge);
	    Py_INCREF(nodeOther);
	    CxEdgeDetach(edge);
	    CxEdgeAttach(edge, node, nodeOther);
	    Py_DECREF(nodeOther);
	    Py_DECREF(edge);
	}
	Py_DECREF(node);
    }

    return rVal;
}

struct CxsTreeCanonizeMinNodeCallback
{
    CxtNodeObject *minNode;
    unsigned taxonNum;
};

static bool
CxpTreeCanonizeNodeCallback(CxtNodeObject *aNode, CxtTreeIteratorStage aStage,
			    void *aContext)
{
    struct CxsTreeCanonizeMinNodeCallback *context
	= (struct CxsTreeCanonizeMinNodeCallback *) aContext;
    unsigned taxonNum;

    // Don't bother looking at nodes twice.
    if (aStage == CxTreeIteratorStagePre)
    {
	if ((context->minNode == NULL)
	    || ((taxonNum = CxNodeTaxonNumGet(aNode)) != CxmTrNodeTaxonNone
		&& (taxonNum < context->taxonNum
		    || context->taxonNum == CxmTrNodeTaxonNone)))
	{
	    context->minNode = aNode;
	    context->taxonNum = taxonNum;
	}
    }

    return false;
}

PyObject *
CxTreeCanonize(CxtTreeObject *self)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;

    CxmXepBegin();
    CxmXepTry
    {
	CxtNodeObject *base;
	CxtRingObject *ring;

	if ((base = CxTreeBaseGet(self)) != NULL)
	{
	    struct CxsTreeCanonizeMinNodeCallback context
		= {NULL, CxmTrNodeTaxonNone};

	    // Set base to be the lowest-numbered taxon.
	    CxTreeIterate(self, CxpTreeCanonizeNodeCallback, NULL, NULL,
			  &context);
	    if (context.minNode != base)
	    {
		base = context.minNode;
		CxTreeBaseSet(self, base);
	    }

	    ring = CxNodeRing(base);
	    if (ring != NULL)
	    {
		CxpTreeCanonize(self, CxRingOther(ring));
	    }
	}

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
