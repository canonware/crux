//==============================================================================
//
// <Copyright = jasone>
// <License>
//
//==============================================================================
//
// Version: <Version = crux>
//
//==============================================================================

#include "../include/_cruxmodule.h"

// Collapse the edge that aRing is part of.  At the end of this function, aRing,
// edge, ringOther, and nodeOther will have been removed from the tree.
// Following is a diagram of how variables are related just before detaching an
// edgeTemp.
//
//                             nodeTemp
//                               |
//                               |
//                             [ring]
//                               |
//                               |
//                           edgeTemp
//                               |
//                               |
//                           ringTemp...
//                          /
//                         /
// aNode--aRing--edge--ringOther--nodeOther
//        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//                  ^^        ^^
//                  || Remove ||
static void
CxpTreeCollapse(CxtNodeObject *aNode, CxtRingObject *aRing)
{
    CxtNodeObject *nodeOther, *nodeTemp;
    CxtEdgeObject *edge, *edgeTemp;
    CxtRingObject *ringOther, *ringTemp;

    edge = CxRingEdge(aRing);
    ringOther = CxRingOther(aRing);
    nodeOther = CxRingNode(ringOther);

    // Hold an extra ref to nodeOther for the duration of rearrangement.  We
    // don't want it to be destroyed until everything has been detached from it.
    Py_INCREF(nodeOther);

    // Iteratively detach edges/nodes from nodeOther, and attach them to node.
    CxmAssert(CxNodeDegree(nodeOther) > 1);
    while ((ringTemp = CxRingNext(ringOther)) != ringOther)
    {
	edgeTemp = CxRingEdge(ringTemp);

	// Get node that is on the other end of ringTemp's edge.
	nodeTemp = CxRingNode(CxRingOther(ringTemp));

	// Detach edgeTemp.  Take care to keep the refcounts above zero.
	Py_INCREF(edgeTemp);
	Py_INCREF(nodeTemp);
	CxEdgeDetach(edgeTemp);

	// Reattach nodeTemp to aNode, using edgeTemp.  Once the attachments
	// have been made, the refcounts can be dropped back down.
	CxEdgeAttach(edgeTemp, aNode, nodeTemp);
	Py_DECREF(nodeTemp);
	Py_DECREF(edgeTemp);
    }

    // Drop nodeOther's refcount back down, now that all of its neighbors have
    // been reattached.
    Py_DECREF(nodeOther);

    // Detach edge from aNode/nodeOther.  edge and nodeOther may become
    // immediately inaccessible due to this call.
    CxEdgeDetach(edge);
}

// Recurse through the tree find edges that can be collapsed, then collapse
// them.  Keep track of the total number of edges that are collapsed.
static bool
CxpTreeCollapseRecurse(CxtTreeObject *self, CxtRingObject *aRing,
		       unsigned *rNcollapsed)
{
    bool rVal;
    CxtNodeObject *node;
    CxtRingObject *ring, *ringOther;
    unsigned degree;

    // Get the degree of the node that this ring is a part of.
    node = CxRingNode(aRing);
    degree = CxNodeDegree(node);

    // Recurse if this is an internal node.
    if (degree > 1)
    {
	// Create an array of ring pointers that will contain pointers to rings
	// that will be collapsed.
	CxtRingObject *ringsCollapse[degree - 1];
	unsigned i, nRingsCollapse = 0;

	// Find edges that should be collapsed and make a list of them.
	for (ring = CxRingNext(aRing); ring != aRing; ring = CxRingNext(ring))
	{
	    ringOther = CxRingOther(ring);
	    if (CxpTreeCollapseRecurse(self, ringOther, rNcollapsed))
	    {
		if (CxEdgeLengthGet(CxRingEdge(ring)) <= 0.0)
		{
		    // Make a note to collapse the edge that ring is part of.
		    // Don't do it until later though, since edge collapsing
		    // will change aNode's ring such that iterative recursion
		    // would be confused about which subtrees had already been
		    // recursed into.
		    ringsCollapse[nRingsCollapse] = ring;
		    nRingsCollapse++;
		}
	    }
	}

	// Collapse edges, now that all collapsable edges have been found.
	for (i = 0; i < nRingsCollapse; i++)
	{
	    CxpTreeCollapse(node, ringsCollapse[i]);
	    (*rNcollapsed)++;
	}

	rVal = true;
    }
    else
    {
	rVal = false;
    }

    return rVal;
}

// Collapse negative- and zero-length edges by creating polytomies.
PyObject *
CxTreeCollapse(CxtTreeObject *self)
{
    CxtNodeObject *base;
    unsigned ncollapsed = 0;

    if ((base = CxTreeBaseGet(self)) != NULL)
    {
	CxtRingObject *ring;

	// Get base's ring.
	if ((ring = CxNodeRing(base)) != NULL)
	{
	    // Recurse.  Collapse negative/zero-length branch only if base isn't
	    // a leaf node, and its neighbor isn't a leaf node.
	    if (CxpTreeCollapseRecurse(self, CxRingOther(ring), &ncollapsed)
		&& CxNodeDegree(base) > 1)
	    {
		// XXX Isn't it okay for base to be a leaf node?
		CxpTreeCollapse(base, ring);
		ncollapsed++;
	    }
	}
    }

    return Py_BuildValue("i", ncollapsed);
}
