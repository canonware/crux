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

typedef struct CxsTreeTbrBisection CxtTreeTbrBisection;
typedef struct CxsTreeTbrData CxtTreeTbrData;

struct CxsTreeTbrBisection
{
    // Bisection edge.
    CxtEdgeObject *edge;

    // Number of TBR neighbors that come before those for this bisection edge.
    unsigned offset;

    // setElm is used to store an edge in the set of edges on each side of a
    // bisection.  This is only incidentally stored in this structure; it's
    // merely convenient to do so.
    CxtEdgeObject *setElm;
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
    unsigned nSetB;
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
		      CxtEdgeObject **rEdgeSet, unsigned *rNEdges)
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

    CxpTreeTbrBEdgeSetGen(self, aRingA, aData->edgeSets, &aData->nSetA);
    CxpTreeTbrBEdgeSetGen(self, aRingB, &aData->edgeSets[aData->nSetA],
			  &aData->nSetA);

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
	    // Record offset.
	    data->bisections[i].offset = offset;

	    // Update offset.
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
		offset += n;
		j++;
	    }
	}
	// Store total number of neighbors in element past the end of the
	// bisection edge list.
	data->bisections[i].offset = offset;

	data->seq = CxTreeSeq(self);
    }

    *rData = data;
    rVal = false;
    RETURN:
    return rVal;
}

bool
CxTreeTbrBEdgeSetsGet(CxtTreeObject *self, CxtEdgeObject *aEdge,
		      CxtEdgeObject ***rSetA, unsigned *rNSetA,
		      CxtEdgeObject ***rSetB, unsigned *rNSetB)
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
    }
    else
    {
	*rSetA = NULL;
    }
    *rNSetA = data->nSetA;

    if (data->nSetB != 0)
    {
	*rSetB = &data->edgeSets[data->nSetA];
    }
    else
    {
	*rSetB = NULL;
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

    if (CxpTreeTbrUpdate(self, &data))
    {
	rVal = true;
	goto RETURN;
    }

//     CxTrTbr(self->tr, aBisect, aReconnectA, aReconnectB);

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

    if (PyArg_ParseTuple(args, "(O!O!O!)",
			 &CxtEdge, &bisect,
			 &CxtEdge, &reconnectA,
			 &CxtEdge, &reconnectB)
	== 0)
    {
	rVal = NULL;
	goto RETURN;
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
CxTreeTbrNneighborsGet(CxtTreeObject *self, unsigned *rNneighbors)
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
	*rNneighbors = 0;
    }
    else
    {
	*rNneighbors = data->bisections[data->nBisections].offset;
    }

    rVal = false;
    RETURN:
    return rVal;
}

PyObject *
CxTreeTbrNneighborsGetPargs(CxtTreeObject *self)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    unsigned nneighbors;

    if (CxTreeTbrNneighborsGet(self, &nneighbors))
    {
	rVal = PyErr_NoMemory();
	goto RETURN;
    }

    rVal = Py_BuildValue("I", nneighbors);
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

    if (CxpTreeTbrUpdate(self, &data))
    {
	rVal = true;
	goto RETURN;
    }

//     CxTrTbrNeighborGet(self->tr, neighbor,
// 		       &bisect, &reconnectA, &reconnectB);

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

    if (PyArg_ParseTuple(args, "I", &neighbor) == 0)
    {
	rVal = NULL;
	goto RETURN;
    }

    if (CxTreeTbrNneighborsGet(self, &nneighbors))
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

    rVal = Py_BuildValue("(O!O!O!)",
			 &CxtEdge, bisect,
			 &CxtEdge, reconnectA,
			 &CxtEdge, reconnectB);
    RETURN:
    return rVal;
}
