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

//#define CxmTreeGCVerbose
#ifdef CxmTreeGCVerbose
#undef Py_INCREF
#define Py_INCREF(op)							\
	fprintf(stderr, "%s:%d:%s(): INCREF(%p) --> %d\n",		\
		__FILE__, __LINE__, __func__, op,			\
		((op)->ob_refcnt) + 1);					\
	(_Py_INC_REFTOTAL  _Py_REF_DEBUG_COMMA				\
	(op)->ob_refcnt++)

#undef Py_DECREF
#define Py_DECREF(op)							\
	fprintf(stderr, "%s:%d:%s(): DECREF(%p) --> %d\n",		\
		__FILE__, __LINE__, __func__, op,			\
		((op)->ob_refcnt) - 1);					\
	if (_Py_DEC_REFTOTAL  _Py_REF_DEBUG_COMMA			\
	   --(op)->ob_refcnt != 0)					\
		_Py_CHECK_REFCNT(op)					\
	else								\
		_Py_Dealloc((PyObject *)(op))
#endif

typedef struct CxsTreeTbrBisection CxtTreeTbrBisection;
typedef struct CxsTreeTbrData CxtTreeTbrData;

struct CxsTreeTbrBisection
{
    // Bisection edge.
    CxtEdgeObject *edge;

    // Number of TBR neighbors that come before those for this bisection edge.
    unsigned offset;
};

struct CxsTreeTbrData
{
    uint64_t seq;

#define CxmTreeTbrInitNBisections 32
    CxtTreeTbrBisection *bisections;
    unsigned maxBisections;
    unsigned nBisections;

    //
    // Bisection edge set data.
    //
    CxtEdgeObject **edgeSets;
    unsigned maxEdgeSets;
    unsigned nSetA;
    CxtRingObject *ringA; // Only valid if nSetA is 0.
    unsigned nSetB;
    CxtRingObject *ringB; // Only valid if nSetB is 0.
};

static void
CxpTreeTbrCleanupFinal(CxtTreeObject *aTree, void *aData, unsigned aInd)
{
    CxtTreeTbrData *data = (CxtTreeTbrData *) aData;

    if (data->bisections != NULL)
    {
	free(data->bisections);
    }

    free(data);
}

static bool
CxpTreeTbrBisectEdgeCallback(CxtEdgeObject *aEdge, CxtTreeIteratorStage aStage,
			     void *aData)
{
    bool rVal;

    if (aStage == CxTreeIteratorStagePre)
    {
	CxtTreeTbrData *data = (CxtTreeTbrData *) aData;

	// Make sure there is enough space to record aEdge.  There must be one
	// extra element in which is stored the total number of TBR neighbors.
	// This is important for the function of bsearch(3).
	if (data->nBisections + 1 < data->maxBisections)
	{
	    // Do nothing.
	}
	else if (data->nBisections == 0)
	{
	    // Allocate for the first time.
	    data->bisections
		= (CxtTreeTbrBisection *)
		malloc(CxmTreeTbrInitNBisections * sizeof(CxtTreeTbrBisection));
	    if (data->bisections == NULL)
	    {
		rVal = true;
		goto RETURN;
	    }

	    data->maxBisections = CxmTreeTbrInitNBisections;
	}
	else
	{
	    CxtTreeTbrBisection *tBisections;

	    // Reallocate.
	    tBisections
		= (CxtTreeTbrBisection *)
		realloc(data->bisections,
			data->maxBisections << 1
			* sizeof(CxtTreeTbrBisection));
	    if (tBisections == NULL)
	    {
		rVal = true;
		goto RETURN;
	    }
	    data->bisections = tBisections;
	    data->maxBisections <<= 1;
	}

	// Record edge.
	data->bisections[data->nBisections].edge = aEdge;
	data->nBisections++;
    }

    rVal = false;
    RETURN:
    return rVal;
}

static void
CxpTreeTbrBEdgeSetGenRecurse(CxtTreeObject *self, CxtRingObject *aRing,
			     CxtEdgeObject **rEdgeSet, unsigned *rNEdges)
{
    CxtRingObject *ring;

    for (ring = CxRingNext(aRing); ring != aRing; ring = CxRingNext(ring))
    {
	// Add edge to set.
	rEdgeSet[*rNEdges] = CxRingEdge(ring);
	(*rNEdges)++;

	// Subtree.
	CxpTreeTbrBEdgeSetGenRecurse(self, CxRingOther(ring), rEdgeSet,
				     rNEdges);
    }
}

CxmpInline void
CxpTreeTbrBEdgeSetGen(CxtTreeObject *self, CxtRingObject *aRing,
		      CxtEdgeObject **rEdgeSet, unsigned *rNEdges,
		      CxtRingObject **rRing)
{
    // Initialize the length of the list before recursing.
    *rNEdges = 0;

    switch (CxNodeDegree(CxRingNode(aRing)))
    {
	case 1:
	{
	    // A subtree that is composed of a single node has no edges.  Add a
	    // single entry to the list, and return the node.
	    rEdgeSet[0] = NULL;
	    (*rNEdges)++;
	    *rRing = aRing;
	    break;
	}
	case 2:
	{
	    // A tree should never have nodes of degree 2.
	    CxmNotReached();
	    break;
	}
	case 3:
	{
	    CxtRingObject *ring;

	    // Take care to add only one of the edges that is connected to the
	    // node, since from the perspective of TBR, the node does not exist.
	    // (A node of degree 2 is a superfluous node.)

	    // First edge.
	    ring = CxRingNext(aRing);
	    rEdgeSet[0] = CxRingEdge(ring);
	    (*rNEdges)++;

	    // First subtree.
	    CxpTreeTbrBEdgeSetGenRecurse(self, CxRingOther(ring), rEdgeSet,
					 rNEdges);

	    // Second subtree.
	    ring = CxRingNext(ring);
	    CxpTreeTbrBEdgeSetGenRecurse(self, CxRingOther(ring), rEdgeSet,
					 rNEdges);

	    break;
	}
	default:
	{
	    CxtRingObject *ring;

	    // Add all edges in the subtree.  Removing the bisection edge still
	    // leaves enough edges attached to the node for the node to have
	    // relevance.

	    for (ring = CxRingNext(aRing);
		 ring != aRing;
		 ring = CxRingNext(ring))
	    {
		// Add edge to set.
		rEdgeSet[*rNEdges] = CxRingEdge(ring);
		(*rNEdges)++;

		// Subtree.
		CxpTreeTbrBEdgeSetGenRecurse(self, CxRingOther(ring), rEdgeSet,
					     rNEdges);
	    }

	    break;
	}
    }
}

static bool
CxpTreeTbrBEdgeSetsGen(CxtTreeObject *self, CxtTreeTbrData *aData,
		       CxtRingObject *aRingA, CxtRingObject *aRingB)
{
    bool rVal;

    // Make sure that there is adequate space for the edge sets.
    if (aData->maxEdgeSets < aData->nBisections)
    {
	if (aData->maxEdgeSets == 0)
	{
	    aData->edgeSets
		= (CxtEdgeObject **) malloc(aData->nBisections
					    * sizeof(CxtEdgeObject *));
	    if (aData->edgeSets == NULL)
	    {
		rVal = true;
		goto RETURN;
	    }
	}
	else
	{
	    CxtEdgeObject **tEdgeSets;

	    tEdgeSets = (CxtEdgeObject **) realloc(aData->edgeSets,
						   aData->nBisections
						   * sizeof(CxtEdgeObject *));
	    if (tEdgeSets == NULL)
	    {
		rVal = true;
		goto RETURN;
	    }

	    aData->edgeSets = tEdgeSets;
	}
	aData->maxEdgeSets = aData->nBisections;
    }

    CxpTreeTbrBEdgeSetGen(self, aRingA, aData->edgeSets, &aData->nSetA,
			  &aData->ringA);
    CxpTreeTbrBEdgeSetGen(self, aRingB, &aData->edgeSets[aData->nSetA],
			  &aData->nSetB, &aData->ringB);

    rVal = false;
    RETURN:
    return rVal;
}

static bool
CxpTreeTbrUpdate(CxtTreeObject *self, CxtTreeTbrData **rData)
{
    bool rVal;
    unsigned treeAuxInd;
    CxtTreeTbrData *data;

    // Get aux indices for TBR data.
    if (CxTreeAuxSearch(self, "TBR", &treeAuxInd))
    {
	// No aux registration.
	data = (CxtTreeTbrData *) malloc(sizeof(CxtTreeTbrData));
	if (data == NULL)
	{
	    rVal = true;
	    goto RETURN;
	}

	data->seq = 0;
	data->bisections = NULL;
	data->maxBisections = 0;
	data->nBisections = 0;

	data->edgeSets = NULL;
	data->maxEdgeSets = 0;
	data->nSetA = 0;
	data->nSetB = 0;

	if (CxTreeAuxRegister(self, "TBR", (void *) data,
			      CxpTreeTbrCleanupFinal, NULL, &treeAuxInd))
	{
	    free(data);
	    rVal = true;
	    goto RETURN;
	}
    }
    else
    {
	data = (CxtTreeTbrData *) CxTreeAuxData(self, treeAuxInd);
    }

    // Update data if the tree has changed, or if the TBR data have never been
    // initialized.
    if (data->seq != CxTreeSeq(self))
    {
	unsigned i, j, offset, n, degreeA, degreeB;
	CxtRingObject *ringA, *ringB;

	// Traverse the tree, and initialize data->bisections along the way.
	if (CxTreeIterate(self, NULL, CxpTreeTbrBisectEdgeCallback, NULL, data))
	{
	    rVal = true;
	    goto RETURN;
	}

	// Iteratively fill in per-edge data.
	for (i = j = offset = 0; i < data->nBisections; i++)
	{
	    // Calculate number of neighbors reachable from this bisection.
	    CxEdgeRingsGet(data->bisections[i].edge, &ringA, &ringB);
	    if (CxpTreeTbrBEdgeSetsGen(self, data, ringA, ringB))
	    {
		rVal = true;
		goto RETURN;
	    }

	    n = (data->nSetA * data->nSetB);
	    // TBR can only be reversed for an edge that is attached to two
	    // nodes of degree 3.
	    degreeA = CxNodeDegree(CxRingNode(ringA));
	    degreeB = CxNodeDegree(CxRingNode(ringB));
	    if ((degreeA == 1 || degreeA == 3)
		&& (degreeB == 1 || degreeB == 3))
	    {
		n--;
	    }

	    if (n != 0)
	    {
		data->bisections[j].edge = data->bisections[i].edge;
		data->bisections[j].offset = offset;
		offset += n;
		j++;
	    }
	}
	// Store total number of neighbors in element past the end of the
	// bisection edge list.
	data->bisections[j].offset = offset;

	// It may be that not all bisections result in neighbors, so the table
	// may not be full.  Keep track of the number of valid elements (not
	// counting the trailing one that stores the total number of TBR
	// neighbors).
	data->nBisections = j;

	data->seq = CxTreeSeq(self);
    }

    *rData = data;
    rVal = false;
    RETURN:
    return rVal;
}

// As part of TBR, extract a node that has only two neighbors.  Take care to
// leave reconnection edges in the tree.  Return NULL, unless there is only one
// node in the subtree; in that case, return the node so that it can be used
// directly during reconnection.
CxmpInline CxtNodeObject *
CxpTreeTbrNodeExtract(CxtTreeObject *self, CxtNodeObject *aNode,
		      CxtEdgeObject *aReconnectA, CxtEdgeObject *aReconnectB,
		      CxtEdgeObject **arTEdges, unsigned *arNTEdges,
		      CxtNodeObject **arTNodes, unsigned *arNTNodes)
{
    CxtNodeObject *rVal;

    switch (CxNodeDegree(aNode))
    {
	case 0:
	{
	    // This node is the only node remaining in the subtree.  It must be
	    // directly reconnected to, so return it.
	    rVal = aNode;
	    break;
	}
	case 1:
	{
	    CxmNotReached();
	}
	case 2:
	{
	    CxtRingObject *ringA, *ringB;
	    CxtRingObject *ringAOther, *ringBOther;
	    CxtEdgeObject *edgeA, *edgeB;
	    CxtNodeObject *nodeA, *nodeB;
	    CxtRingObject *tRingA, *tRingB;

	    // Get all variables that are necessary for careful extraction of
	    // aNode, and proper rematching of rings with nodes.  The
	    // rematching is critical to the maintenance of the character state
	    // sets in leaf nodes (which node[AB] may or may not be).
	    ringA = CxNodeRing(aNode);
	    edgeA = CxRingEdge(ringA);
	    ringAOther = CxRingOther(ringA);
	    nodeA = CxRingNode(ringAOther);

	    ringB = CxRingNext(ringA);
	    edgeB = CxRingEdge(ringB);
	    ringBOther = CxRingOther(ringB);
	    nodeB = CxRingNode(ringBOther);

	    // Detach.
	    Py_INCREF(aNode); // +1
	    Py_INCREF(edgeA); // +1
	    Py_INCREF(nodeA); // +1
	    CxEdgeDetach(edgeA);
	    Py_INCREF(nodeB); // +1
	    Py_INCREF(edgeB); // +1
	    CxEdgeDetach(edgeB);

	    // Store aNode as a spare.
	    arTNodes[*arNTNodes] = aNode;
	    (*arNTNodes)++;

	    // Be careful to preserve reconnection edges, which either edgeA or
	    // edgeB may be.
	    if (edgeB != aReconnectA && edgeB != aReconnectB)
	    {
		// Use edgeA when splicing node[AB] back together.

		// Swap data in ringA and ringBOther.
		CxRingAuxSwap(ringA, ringBOther);

		// Attach node[AB].  Take care to keep the proper ends of edgeA
		// associated with node[AB].
		CxEdgeRingsGet(edgeA, &tRingA, &tRingB);
		if (ringAOther == tRingA)
		{
		    CxEdgeAttach(edgeA, nodeA, nodeB);
		}
		else
		{
		    CxEdgeAttach(edgeA, nodeB, nodeA);
		}
		Py_DECREF(edgeA); // -1

		// Store edgeB as a spare.
		arTEdges[*arNTEdges] = edgeB;
		(*arNTEdges)++;
	    }
	    else
	    {
		// Use edgeB when splicing node[AB] back together.
		CxmAssert(edgeA != aReconnectA && edgeA != aReconnectB);

		// Swap data in ringB and ringAOther.
		CxRingAuxSwap(ringB, ringAOther);

		// Attach node[AB].  Take care to keep the proper ends of
		// edgeB associated with node[AB].
		CxEdgeRingsGet(edgeB, &tRingA, &tRingB);
		if (ringB == tRingA)
		{
		    CxEdgeAttach(edgeB, nodeA, nodeB);
		}
		else
		{
		    CxEdgeAttach(edgeB, nodeB, nodeA);
		}
		Py_DECREF(edgeB); // -1

		// Store edgeA as a spare.
		arTEdges[*arNTEdges] = edgeA;
		(*arNTEdges)++;
	    }

	    Py_DECREF(nodeB); // -1
	    Py_DECREF(nodeA); // -1

	    rVal = NULL;
	    break;
	}
	default:
	{
	    // Do nothing, since this node has enough neighbors to remain
	    // relevant (3 or more).
	    rVal = NULL;
	}
    }

    return rVal;
}

// Splice a node into the middle of aEdge, and return the node.
CxmpInline CxtNodeObject *
CxpTreeTbrNodeSplice(CxtTreeObject *self, CxtEdgeObject *aEdge,
		     CxtEdgeObject **arTEdges, unsigned *arNTEdges,
		     CxtNodeObject **arTNodes, unsigned *arNTNodes)
{
    CxtNodeObject *rVal;
    CxtNodeObject *nodeA, *nodeB;
    CxtRingObject *ringA, *ringB;
    CxtEdgeObject *spliceEdge;
    CxtRingObject *spliceRingA, *spliceRingB;

    // Get all variables that are necessary for careful splicing of a node into
    // aEdge, and proper rematching of rings with nodes.  The rematching is
    // critical to the MP code (maintenance of the character state sets in leaf
    // nodes, which node[AB] may or may not be).
    CxEdgeRingsGet(aEdge, &ringA, &ringB);
    nodeA = CxRingNode(ringA);
    nodeB = CxRingNode(ringB);

    // Get an edge.
    CxmAssert(*arNTEdges > 0);
    (*arNTEdges)--;
    spliceEdge = arTEdges[*arNTEdges];
    CxEdgeRingsGet(spliceEdge, &spliceRingA, &spliceRingB);

    // Get a node.
    CxmAssert(*arNTNodes > 0);
    (*arNTNodes)--;
    rVal = arTNodes[*arNTNodes];

    // Detach.
    Py_INCREF(nodeA); // +1
    Py_INCREF(nodeB); // +1
    CxEdgeDetach(aEdge);

    // Swap data in ringB and spliceRingA.
    CxRingAuxSwap(ringB, spliceRingA);

    // Reattach.
    CxEdgeAttach(aEdge, nodeA, rVal);
    CxEdgeAttach(spliceEdge, nodeB, rVal);
    Py_DECREF(nodeB); // -1
    Py_DECREF(nodeA); // -1
    Py_DECREF(spliceEdge); // -1

    return rVal;
}

// Comparison function that is passed to bsearch(3).
static int
CxpTreeTbrBisectionCompare(const void *aKey, const void *aVal)
{
    int rVal;
    const CxtTreeTbrBisection *key = (const CxtTreeTbrBisection *) aKey;
    const CxtTreeTbrBisection *val = (const CxtTreeTbrBisection *) aVal;

    if (key->offset < val->offset)
    {
	rVal = -1;
    }
    else if (key->offset < (&val[1])->offset)
    {
	rVal = 0;
    }
    else
    {
	rVal = 1;
    }

    return rVal;
}

bool
CxTreeTbrBEdgeSetsGet(CxtTreeObject *self, CxtEdgeObject *aEdge,
		      CxtEdgeObject ***rSetA, unsigned *rNSetA,
		      CxtRingObject **rRingA,
		      CxtEdgeObject ***rSetB, unsigned *rNSetB,
		      CxtRingObject **rRingB)
{
    bool rVal;
    CxtTreeTbrData *data;
    CxtRingObject *ringA, *ringB;

    CxEdgeRingsGet(aEdge, &ringA, &ringB);
    if (CxpTreeTbrUpdate(self, &data)
	|| CxpTreeTbrBEdgeSetsGen(self, data, ringA, ringB))
    {
	rVal = true;
	goto RETURN;
    }

    if (data->nSetA != 0)
    {
	*rSetA = data->edgeSets;
	*rRingA = NULL;
    }
    else
    {
	*rSetA = NULL;
	*rRingA = data->ringA;
    }
    *rNSetA = data->nSetA;

    if (data->nSetB != 0)
    {
	*rSetB = &data->edgeSets[data->nSetA];
	*rRingB = NULL;
    }
    else
    {
	*rSetB = NULL;
	*rRingB = data->ringB;
    }
    *rNSetB = data->nSetB;

    rVal = false;
    RETURN:
    return rVal;
}

bool
CxTreeTbr(CxtTreeObject *self, CxtEdgeObject *aBisect,
	  CxtEdgeObject *aReconnectA, CxtEdgeObject *aReconnectB)
{
    bool rVal;
    CxtTreeTbrData *data;
    CxtRingObject *ringA, *ringB;
    CxtNodeObject *nodeA, *nodeB, *nodes[4] = {NULL, NULL, NULL, NULL};
    CxtEdgeObject *tEdges[2];
    unsigned nTEdges = 0;
    CxtNodeObject *tNodes[2];
    unsigned nTNodes = 0;

    if (CxpTreeTbrUpdate(self, &data))
    {
	rVal = true;
	goto RETURN;
    }

    // Get the nodes to either side of the edge where the bisection will be
    // done.
    CxEdgeRingsGet(aBisect, &ringA, &ringB);
    nodeA = CxRingNode(ringA);
    nodeB = CxRingNode(ringB);

    // Determine whether additional edges and nodes will be needed for
    // reconnection.  If so, allocate them here, so that if there is an OOM
    // error, recovery is simpler.
    if (CxNodeDegree(nodeA) > 3)
    {
	tEdges[nTEdges] = CxEdgeNew(self); // =1
	if (tEdges[nTEdges] == NULL)
	{
	    rVal = true;
	    goto RETURN;
	}
	nTEdges++;

	tNodes[nTNodes] = CxNodeNew(self); // =1
	if (tNodes[nTNodes] == NULL)
	{
	    Py_DECREF(tEdges[0]); // -1
	    rVal = true;
	    goto RETURN;
	}
	nTNodes++;
    }

    if (CxNodeDegree(nodeB) > 3)
    {
	tEdges[nTEdges] = CxEdgeNew(self); // =1
	if (tEdges[nTEdges] == NULL)
	{
	    Py_DECREF(tNodes[0]); // -1
	    Py_DECREF(tEdges[0]); // -1
	    rVal = true;
	    goto RETURN;
	}
	nTEdges++;

	tNodes[nTNodes] = CxNodeNew(self); // =1
	if (tNodes[nTNodes] == NULL)
	{
	    Py_DECREF(tEdges[1]); // -1
	    Py_DECREF(tNodes[0]); // -1
	    Py_DECREF(tEdges[0]); // -1
	    rVal = true;
	    goto RETURN;
	}
	nTNodes++;
    }

    // Bisect.  aBisect will be used below for reconnection.
    Py_INCREF(nodeA); // +1
    Py_INCREF(nodeB); // +1
    Py_INCREF(aBisect); // +1
    CxEdgeDetach(aBisect);

    // For node[AB], extract the node if it has only two neighbors.
    //
    // nodes[0..1] are NULL, unless they refer to the only node in a subtree.
    nodes[0] = CxpTreeTbrNodeExtract(self, nodeA, aReconnectA, aReconnectB,
				     tEdges, &nTEdges, tNodes, &nTNodes);
    // +1? nodeA === tNodes[nTNodes - 1]
    // +1? tEdges[nTEdges - 1]
    CxmAssert(nTEdges <= 2);
    CxmAssert(nTNodes <= 2);

    nodes[1] = CxpTreeTbrNodeExtract(self, nodeB, aReconnectA, aReconnectB,
				     tEdges, &nTEdges, tNodes, &nTNodes);
    // +1? nodeB === tNodes[nTNodes - 1]
    // +1? tEdges[nTEdges - 1]
    CxmAssert(nTEdges <= 2);
    CxmAssert(nTNodes <= 2);

    // For each reconnection edge, splice a node into the edge (if the subtree
    // has more than one node).
    //
    // nodes[2..3] are set to NULL if no reconnection edge is specified.
    if (aReconnectA != NULL)
    {
	nodes[2] = CxpTreeTbrNodeSplice(self, aReconnectA,
					tEdges, &nTEdges, tNodes, &nTNodes);
	// +1? nodes[2] === -1? tNodes[nTNodes]
	// -1 tEdges[nTEdges]
    }
    else
    {
	nodes[2] = NULL;
    }

    if (aReconnectB != NULL)
    {
	nodes[3] = CxpTreeTbrNodeSplice(self, aReconnectB,
					tEdges, &nTEdges, tNodes, &nTNodes);
	// +1? nodes[3] === -1? tNodes[nTNodes]
	// -1 tEdges[nTEdges]
    }
    else
    {
	nodes[3] = NULL;
    }

    CxmAssert(nTEdges == 0);
    CxmAssert(nTNodes == 0);

    // If either subtree has only a single node, special care must be taken
    // during reconnection to re-associate the proper end of the bisection edge
    // with the single node.  This is primarily because for MP character state
    // information for leaf nodes is actually stored in the rings that attach
    // them to the tree, and breaking this association would require
    // re-initializing the character state set vectors.
    if (nodes[0] != NULL)
    {
	// nodes[0] (same as nodeA) is a single node.
	CxmAssert(nodes[0] == nodeA);
	CxmAssert(nodes[1] == NULL);

	if (nodes[2] == NULL)
	{
	    nodes[2] = nodes[3];
	    nodes[3] = NULL;
	}
	CxEdgeAttach(aBisect, nodes[0], nodes[2]);
    }
    else if (nodes[1] != NULL)
    {
	// nodes[1] (same as nodeB) is a single node.
	CxmAssert(nodes[1] == nodeB);

	if (nodes[2] == NULL)
	{
	    nodes[2] = nodes[3];
	    nodes[3] = NULL;
	}
	CxEdgeAttach(aBisect, nodes[2], nodes[1]);
    }
    else
    {
	// Bisection was not done adjacent to a leaf node.  Attach the two
	// spliced-in nodes.
	CxEdgeAttach(aBisect, nodes[2], nodes[3]);
    }

    // Clean up references.
    if (nodes[3] != NULL)
    {
	Py_DECREF(nodes[3]); // -1
    }
    if (nodes[2] != NULL)
    {
	Py_XDECREF(nodes[2]); // -1
    }

    Py_DECREF(aBisect); // -1
    Py_DECREF(nodeB); // -1
    Py_DECREF(nodeA); // -1

    rVal = false;
    RETURN:
    return rVal;
}

PyObject *
CxTreeTbrPargs(CxtTreeObject *self, PyObject *args)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    CxtEdgeObject *bisect, *reconnectA, *reconnectB;

    if (PyArg_ParseTuple(args, "(O!OO)",
			 &CxtEdge, &bisect,
			 &reconnectA,
			 &reconnectB)
	== 0)
    {
	rVal = NULL;
	goto RETURN;
    }

    // Convert None reconnection edges to NULL, if necessary.  Otherwise, make
    // sure that the arguments are edge objects.
    if (reconnectA != (CxtEdgeObject *) Py_None)
    {
	if (PyObject_TypeCheck(reconnectA, &CxtEdge) == 0)
	{
	    CxError(CxgTreeTypeError, "reconnectA: Edge or None expected");
	    rVal = NULL;
	    goto RETURN;
	}
    }
    else
    {
	reconnectA = NULL;
    }

    if (reconnectB != (CxtEdgeObject *) Py_None)
    {
	if (PyObject_TypeCheck(reconnectB, &CxtEdge) == 0)
	{
	    CxError(CxgTreeTypeError, "reconnectB: Edge or None expected");
	    rVal = NULL;
	    goto RETURN;
	}
    }
    else
    {
	reconnectB = NULL;
    }

    if (CxTreeTbr(self, bisect, reconnectA, reconnectB))
    {
	rVal = PyErr_NoMemory();
	goto RETURN;
    }

    Py_INCREF(Py_None);
    rVal = Py_None;
    RETURN:
    return rVal;
}

bool
CxTreeTbrNEdgesGet(CxtTreeObject *self, unsigned *rNEdges)
{
    bool rVal;
    CxtTreeTbrData *data;

    if (CxpTreeTbrUpdate(self, &data))
    {
	rVal = true;
	goto RETURN;
    }

    *rNEdges = data->nBisections;

    rVal = false;
    RETURN:
    return rVal;
}

bool
CxTreeTbrEdgeGet(CxtTreeObject *self, unsigned aEdge, CxtEdgeObject **rEdge)
{
    bool rVal;
    CxtTreeTbrData *data;

    if (CxpTreeTbrUpdate(self, &data))
    {
	rVal = true;
	goto RETURN;
    }

    CxmAssert(aEdge < data->nBisections);
    *rEdge = data->bisections[aEdge].edge;

    rVal = false;
    RETURN:
    return rVal;
}

bool
CxTreeTbrEdgeOffset(CxtTreeObject *self, unsigned aEdge, unsigned *rOffset)
{
    bool rVal;
    CxtTreeTbrData *data;

    if (CxpTreeTbrUpdate(self, &data))
    {
	rVal = true;
	goto RETURN;
    }

    CxmAssert(aEdge <= data->nBisections);
    *rOffset = data->bisections[aEdge].offset;

    rVal = false;
    RETURN:
    return rVal;
}

bool
CxTreeTbrNNeighborsGet(CxtTreeObject *self, unsigned *rNNeighbors)
{
    bool rVal;
    CxtTreeTbrData *data;

    if (CxpTreeTbrUpdate(self, &data))
    {
	rVal = true;
	goto RETURN;
    }

    if (data->nBisections == 0)
    {
	*rNNeighbors = 0;
    }
    else
    {
	*rNNeighbors = data->bisections[data->nBisections].offset;
    }

    rVal = false;
    RETURN:
    return rVal;
}

PyObject *
CxTreeTbrNNeighborsGetPargs(CxtTreeObject *self)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    unsigned nneighbors;

    if (CxTreeTbrNNeighborsGet(self, &nneighbors))
    {
	rVal = PyErr_NoMemory();
	goto RETURN;
    }

    rVal = Py_BuildValue("i", nneighbors);
    RETURN:
    return rVal;
}

bool
CxTreeTbrNeighborGet(CxtTreeObject *self, unsigned aNeighbor,
		     CxtEdgeObject **rBisect, CxtEdgeObject **rReconnectA,
		     CxtEdgeObject **rReconnectB)
{
    bool rVal;
    CxtTreeTbrData *data;
    CxtTreeTbrBisection key, *bisection;
    CxtRingObject *ringA, *ringB;
    unsigned rem, degreeA, degreeB;

    if (CxpTreeTbrUpdate(self, &data))
    {
	rVal = true;
	goto RETURN;
    }

    // Get the bisection edge.
    CxmAssert(aNeighbor < data->bisections[data->nBisections].offset);
    key.offset = aNeighbor;
    bisection = bsearch(&key, data->bisections, data->nBisections,
			sizeof(CxtTreeTbrBisection),
			CxpTreeTbrBisectionCompare);
    CxmCheckPtr(bisection);
    *rBisect = bisection->edge;

    // Get the reconnection edges.  The indices for a and b are mapped onto the
    // edges of each subtree, in a peculiar fashion.  The actual ordering isn't
    // important, as long as it is consistent (repeatable) and correct (all TBR
    // neighbors are enumerated).  The edge index mapping can be summarized as
    // follows:
    //
    // 1) Start with a full tree.
    //
    // 2) (Pretend to) bisect the tree at the appropriate edge.
    //
    // 3) For each subtree, do a recursive in-order traversal and build a list
    //    of edges, not including one of the edges adjacent to the bisection (in
    //    the case of an internal bifurcating node adjacent to the bisection).
    //
    // 4) Use a nested loop to iterate over neighbors, where each iteration is a
    //    combination of edges in the two subtrees.  The first combination is
    //    skipped, if it would reverse the bisection (if the neighboring nodes
    //    are both of degree 1 or 3).

    // Generate the edge lists.
    CxEdgeRingsGet(bisection->edge, &ringA, &ringB);
    if (CxpTreeTbrBEdgeSetsGen(self, data, ringA, ringB))
    {
	rVal = true;
	goto RETURN;
    }

    // Calculate the offset of the neighbor from the beginning of this edge's
    // reconnection combination enumeration.
    rem = aNeighbor - bisection->offset;

    // Avoid the first combination, if it would reverse the bisection.
    degreeA = CxNodeDegree(CxRingNode(ringA));
    degreeB = CxNodeDegree(CxRingNode(ringB));
    if ((degreeA == 1 || degreeA == 3) && (degreeB == 1 || degreeB == 3))
    {
	rem++;
    }

    *rReconnectA = data->edgeSets[rem / data->nSetB];
    *rReconnectB = data->edgeSets[data->nSetA + (rem % data->nSetB)];
    CxmAssert(*rReconnectA != NULL || *rReconnectB != NULL);

    rVal = false;
    RETURN:
    return rVal;
}

PyObject *
CxTreeTbrNeighborGetPargs(CxtTreeObject *self, PyObject *args)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    unsigned nneighbors, neighbor;
    CxtEdgeObject *bisect, *reconnectA, *reconnectB;

    if (PyArg_ParseTuple(args, "i", &neighbor) == 0)
    {
	rVal = NULL;
	goto RETURN;
    }

    if (CxTreeTbrNNeighborsGet(self, &nneighbors))
    {
	rVal = PyErr_NoMemory();
	goto RETURN;
    }

    if (neighbor >= nneighbors)
    {
	CxError(CxgTreeValueError,
		"neighbor: %u is out of range [0..%u]",
		neighbor, nneighbors);
	rVal = NULL;
	goto RETURN;
    }

    if (CxTreeTbrNeighborGet(self, neighbor, &bisect, &reconnectA, &reconnectB))
    {
	rVal = PyErr_NoMemory();
	goto RETURN;
    }

    // reconnect[AB] may be NULL, which must be translated to None.
    if (reconnectA == NULL)
    {
	Py_INCREF(Py_None);
	reconnectA = (CxtEdgeObject *) Py_None;
    }
    else if (reconnectB == NULL)
    {
	Py_INCREF(Py_None);
	reconnectB = (CxtEdgeObject *) Py_None;
    }

    rVal = Py_BuildValue("(OOO)",
			 (PyObject *) bisect,
			 (PyObject *) reconnectA,
			 (PyObject *) reconnectB);
    RETURN:
    return rVal;
}
