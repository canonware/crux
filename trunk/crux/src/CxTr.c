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
//
// Tr implements multifurcating phylogenetic trees and various operations on
// them.  Nodes can be manipulated via the CxTrNode*() APIs, and edges can be
// manipulated via the CxTrEdge*() APIs.
//
// Nodes and edges are allocated from per-tr internal arrays.  This allows array
// indices to be used for links between nodes and edges.
//
// Each edge has two ring objects associated with it, one for each end of the
// edge.  These ring objects are also allocated from an internal array that has
// a direct correspondence to the edge array.  This provides for constant-time
// conversion between edges and rings.
//
// Each node optionally references a ring, which contains all the edges that
// connect the node to other nodes.  Each ring object in the ring refers to the
// node, which allows constant-time conversion from rings to nodes.
//
//==============================================================================

#define CxmTr_c
#include "../include/_cruxmodule.h"

//==============================================================================

// CxTrRing.

CxmpInline void
CxpTrRingInit(CxtTr *aTr, CxtTrRing aRing)
{
    CxtTrr *trr;

    CxmAssert((aRing >> 1) < aTr->ntres);

    trr = &aTr->trrs[aRing];

    CxmQriNew(aTr->trrs, aRing, link);
    trr->aux = NULL;
    trr->node = CxmTrNodeNone;
    trr->ps = NULL;
}

void *
CxTrRingAuxGet(CxtTr *aTr, CxtTrRing aRing)
{
    CxmDassert(CxpTrRingValidate(aTr, aRing));

    return aTr->trrs[aRing].aux;
}

void
CxTrRingAuxSet(CxtTr *aTr, CxtTrRing aRing, void *aAux)
{
    CxmDassert(CxpTrRingValidate(aTr, aRing));

    aTr->trrs[aRing].aux = aAux;
}

//==============================================================================

// Validation functions.

#ifdef CxmDebug
// Validate a ring.
bool
CxpTrRingValidate(CxtTr *aTr, CxtTrRing aRing)
{
    CxtTrr *trr;

    CxmCheckPtr(aTr);
    CxmAssert(aTr->magic == CxmTrMagic);
    CxmAssert((aRing >> 1) < aTr->ntres);

    trr = &aTr->trrs[aRing];

    CxmAssert(trr->node < aTr->ntrns || trr->node == CxmTrNodeNone);

    return true;
}

// Validate an edge.
bool
CxpTrEdgeValidate(CxtTr *aTr, CxtTrEdge aEdge)
{
    CxtTre *tre;

    CxmCheckPtr(aTr);
    CxmAssert(aTr->magic == CxmTrMagic);
    CxmAssert(aEdge < aTr->ntres);

    tre = &aTr->tres[aEdge];

    CxmAssert(tre->magic == CxmTreMagic);

    CxpTrRingValidate(aTr, (aEdge << 1));
    CxpTrRingValidate(aTr, (aEdge << 1) + 1);

    return true;
}

// Validate a node.
bool
CxpTrNodeValidate(CxtTr *aTr, CxtTrNode aNode)
{
    CxtTrn *trn;
    CxtTrRing ring;
    uint32_t nneighbors;

    CxmCheckPtr(aTr);
    CxmAssert(aTr->magic == CxmTrMagic);
    CxmAssert(aNode < aTr->ntrns);

    trn = &aTr->trns[aNode];

    CxmAssert(trn->magic == CxmTrnMagic);

    nneighbors = 0;
    CxmQliForeach(ring, &trn->rings, aTr->trrs, link)
	{
	    // Validate edge.
	    CxpTrEdgeValidate(aTr, CxTrRingEdgeGet(aTr, ring));

	    nneighbors++;
	}

    if (trn->taxonNum != CxmTrNodeTaxonNone)
    {
	// Only leaf nodes can have taxon numbers.  Leaf nodes have at most 1
	// neighbor.
	CxmAssert(nneighbors <= 1);
    }

    return true;
}
#endif

//==============================================================================

// CxTrPs.

CxmpInline void
CxpTrPsDelete(CxtTr *aTr, CxtTrPs *aPs)
{
    if (aPs->chars != NULL)
    {
	CxmFree(aPs->achars);
    }

    CxmFree(aPs);
}

//==============================================================================

// trEdge.

CxmpInline void
CxpTrEdgeInit(CxtTr *aTr, CxtTrEdge aEdge)
{
    CxtTre *tre;

    tre = &aTr->tres[aEdge];

    tre->u.aux = NULL;
    CxpTrRingInit(aTr, (aEdge << 1));
    CxpTrRingInit(aTr, (aEdge << 1) + 1);
    tre->length = 0.0;
    tre->ps = NULL;

#ifdef CxmDebug
    tre->magic = CxmTreMagic;
#endif
}

CxmpInline CxtTrEdge
CxpTrEdgeAlloc(CxtTr *aTr)
{
    CxtTrEdge rVal;

    if (aTr->sparetres == CxmTrEdgeNone)
    {
	uint32_t i, nspares;

	if (aTr->tres == NULL)
	{
	    aTr->tres = (CxtTre *) CxmMalloc(sizeof(CxtTre));
	    CxmAssert(aTr->trrs == NULL);
	    aTr->trrs = (CxtTrr *) CxmMalloc(sizeof(CxtTrr) * 2);
	    nspares = 1;
	    aTr->ntres = 1;
	}
	else
	{
	    aTr->tres = (CxtTre *) CxmRealloc(aTr->tres,
					      sizeof(CxtTre)
					      * aTr->ntres * 2);
	    CxmCheckPtr(aTr->trrs);
	    aTr->trrs = (CxtTrr *) CxmRealloc(aTr->trrs,
					      sizeof(CxtTrr)
					      * aTr->ntres * 4);
	    nspares = aTr->ntres;
	    aTr->ntres *= 2;
	}

	// Initialize last spare.
	aTr->sparetres = aTr->ntres - 1;
#ifdef CxmDebug
	aTr->tres[aTr->sparetres].magic = 0xa5a5a5a5;
#endif
	aTr->tres[aTr->sparetres].u.link = CxmTrEdgeNone;

	// Insert other spares into spares stack.
	for (i = 1; i < nspares; i++)
	{
	    aTr->sparetres--;
#ifdef CxmDebug
	    aTr->tres[aTr->sparetres].magic = 0xa5a5a5a5;
#endif
	    aTr->tres[aTr->sparetres].u.link = aTr->sparetres + 1;
	}
    }

    // Remove a spare from the spares stack.
    rVal = aTr->sparetres;
    aTr->sparetres = aTr->tres[rVal].u.link;

    // Initialize rVal.
    CxpTrEdgeInit(aTr, rVal);

    return rVal;
}

CxmpInline void
CxpTrEdgeDealloc(CxtTr *aTr, CxtTrEdge aEdge)
{
    CxtTre *tre;
    CxtTrr *trr;

    tre = &aTr->tres[aEdge];
    if (tre->ps != NULL)
    {
	CxpTrPsDelete(aTr, tre->ps);
    }
#ifdef CxmDebug
    memset(tre, 0x5a, sizeof(CxtTre));
#endif

    trr = &aTr->trrs[aEdge << 1];
    if (trr->ps != NULL)
    {
	CxpTrPsDelete(aTr, trr->ps);
    }
#ifdef CxmDebug
    memset(trr, 0x5a, sizeof(CxtTrr));
#endif

    trr = &aTr->trrs[(aEdge << 1) + 1];
    if (trr->ps != NULL)
    {
	CxpTrPsDelete(aTr, trr->ps);
    }
#ifdef CxmDebug
    memset(trr, 0x5a, sizeof(CxtTrr));
#endif

    aTr->tres[aEdge].u.link = aTr->sparetres;
    aTr->sparetres = aEdge;
}

CxtTrEdge
CxTrEdgeNew(CxtTr *aTr)
{
    return CxpTrEdgeAlloc(aTr);
}

void
CxTrEdgeDelete(CxtTr *aTr, CxtTrEdge aEdge)
{
    CxmDassert(CxpTrEdgeValidate(aTr, aEdge));
    CxmAssert(CxTrRingNodeGet(aTr, CxTrEdgeRingGet(aTr, aEdge, 0))
	      == CxmTrNodeNone);
    CxmAssert(CxTrRingNodeGet(aTr, CxTrEdgeRingGet(aTr, aEdge, 1))
	      == CxmTrNodeNone);

    CxpTrEdgeDealloc(aTr, aEdge);
}

double
CxTrEdgeLengthGet(CxtTr *aTr, CxtTrEdge aEdge)
{
    CxmDassert(CxpTrEdgeValidate(aTr, aEdge));

    return aTr->tres[aEdge].length;
}

void
CxTrEdgeLengthSet(CxtTr *aTr, CxtTrEdge aEdge, double aLength)
{
    CxmDassert(CxpTrEdgeValidate(aTr, aEdge));

    aTr->tres[aEdge].length = aLength;
}

void *
CxTrEdgeAuxGet(CxtTr *aTr, CxtTrEdge aEdge)
{
    CxmDassert(CxpTrEdgeValidate(aTr, aEdge));

    return aTr->tres[aEdge].u.aux;
}

void
CxTrEdgeAuxSet(CxtTr *aTr, CxtTrEdge aEdge, void *aAux)
{
    CxmDassert(CxpTrEdgeValidate(aTr, aEdge));

    aTr->tres[aEdge].u.aux = aAux;
}

void
CxTrEdgeAttach(CxtTr *aTr, CxtTrEdge aEdge, CxtTrNode aNodeA,
	       CxtTrNode aNodeB)
{
    CxtTrn *trn;
    CxtTrRing ring;

    CxmDassert(CxpTrEdgeValidate(aTr, aEdge));
    CxmDassert(CxTrRingNodeGet(aTr, CxTrEdgeRingGet(aTr, aEdge, 0))
	       == CxmTrNodeNone);
    CxmDassert(CxTrRingNodeGet(aTr, CxTrEdgeRingGet(aTr, aEdge, 1))
	       == CxmTrNodeNone);
    CxmDassert(CxpTrNodeValidate(aTr, aNodeA));
    CxmDassert(CxpTrNodeValidate(aTr, aNodeB));
    CxmAssert(aNodeA != aNodeB);
    CxmAssert(CxTrNodeDistance(aTr, aNodeA, aNodeB) == 0);

    // First end.
    ring = CxTrEdgeRingGet(aTr, aEdge, 0);
    trn = &aTr->trns[aNodeA];
    CxmQliHeadInsert(&trn->rings, aTr->trrs, ring, link);
    aTr->trrs[ring].node = aNodeA;

    // Second end.
    ring = CxTrEdgeRingGet(aTr, aEdge, 1);
    trn = &aTr->trns[aNodeB];
    CxmQliHeadInsert(&trn->rings, aTr->trrs, ring, link);
    aTr->trrs[ring].node = aNodeB;

    // Mark tree as modified.
    aTr->modified = true;

    CxmDassert(CxpTrEdgeValidate(aTr, aEdge));
    CxmDassert(CxpTrNodeValidate(aTr, aNodeA));
    CxmDassert(CxpTrNodeValidate(aTr, aNodeB));
    CxmAssert(CxTrNodeDistance(aTr, aNodeA, aNodeB) == 1);
}

void
CxTrEdgeDetach(CxtTr *aTr, CxtTrEdge aEdge)
{
    CxtTrRing ring;
    CxtTrn *trn;

    CxmDassert(CxpTrEdgeValidate(aTr, aEdge));
    CxmDassert(CxTrRingNodeGet(aTr, CxTrEdgeRingGet(aTr, aEdge, 0))
	       != CxmTrNodeNone);
    CxmDassert(CxTrRingNodeGet(aTr, CxTrEdgeRingGet(aTr, aEdge, 1))
	       != CxmTrNodeNone);
    CxmAssert(CxTrNodeDistance(aTr,
			       CxTrRingNodeGet(aTr, CxTrEdgeRingGet(aTr, aEdge,
								    0)),
			       CxTrRingNodeGet(aTr, CxTrEdgeRingGet(aTr, aEdge,
								    1))) == 1);

    // Detach from neighboring nodes.  Use CxmQliRemove() to make sure that the
    // nodes still point to their rings.
    ring = CxTrEdgeRingGet(aTr, aEdge, 0);
    trn = &aTr->trns[CxTrRingNodeGet(aTr, ring)];
    CxmQliRemove(&trn->rings, aTr->trrs, ring, link);
    aTr->trrs[ring].node = CxmTrNodeNone;

    ring = CxTrEdgeRingGet(aTr, aEdge, 1);
    trn = &aTr->trns[CxTrRingNodeGet(aTr, ring)];
    CxmQliRemove(&trn->rings, aTr->trrs, ring, link);
    aTr->trrs[ring].node = CxmTrNodeNone;

    // Mark tree as modified.
    aTr->modified = true;

    CxmDassert(CxpTrEdgeValidate(aTr, aEdge));
}

//==============================================================================

// CxTrNode.

CxmpInline void
CxpTrNodeInit(CxtTr *aTr, CxtTrNode aNode)
{
    CxtTrn *trn;

    trn = &aTr->trns[aNode];

    trn->u.aux = NULL;
    trn->taxonNum = CxmTrNodeTaxonNone;
    CxmQliNew(&trn->rings);

#ifdef CxmDebug
    trn->magic = CxmTrnMagic;
#endif
}

CxmpInline CxtTrNode
CxpTrNodeAlloc(CxtTr *aTr)
{
    CxtTrNode rVal;

    if (aTr->sparetrns == CxmTrNodeNone)
    {
	uint32_t i, nspares;

	// Allocate spares.
	if (aTr->trns == NULL)
	{
	    aTr->trns
		= (CxtTrn *) CxmMalloc(sizeof(CxtTrn));
	    nspares = 1;
	    aTr->ntrns = 1;
	}
	else
	{
	    aTr->trns = (CxtTrn *) CxmRealloc(aTr->trns,
					      sizeof(CxtTrn)
					      * aTr->ntrns * 2);
	    nspares = aTr->ntrns;
	    aTr->ntrns *= 2;
	}

	// Initialize last spare.
	aTr->sparetrns = aTr->ntrns - 1;
#ifdef CxmDebug
	aTr->trns[aTr->sparetrns].magic = 0xa5a5a5a5;
#endif
	aTr->trns[aTr->sparetrns].u.link = CxmTrNodeNone;

	// Insert other spares into spares stack.
	for (i = 1; i < nspares; i++)
	{
	    aTr->sparetrns--;
#ifdef CxmDebug
	    aTr->trns[aTr->sparetrns].magic = 0xa5a5a5a5;
#endif
	    aTr->trns[aTr->sparetrns].u.link = aTr->sparetrns + 1;
	}
    }

    // Remove a spare from the spares stack.
    rVal = aTr->sparetrns;
    aTr->sparetrns = aTr->trns[rVal].u.link;

    // Initialize rVal.
    CxpTrNodeInit(aTr, rVal);

    return rVal;
}

CxmpInline void
CxpTrNodeDealloc(CxtTr *aTr, CxtTrNode aNode)
{
#ifdef CxmDebug
    CxtTrn *trn;

    trn = &aTr->trns[aNode];

    memset(trn, 0x5a, sizeof(CxtTrn));
#endif

    aTr->trns[aNode].u.link = aTr->sparetrns;
    aTr->sparetrns = aNode;
}

// Calculate the number of edges connected to the node that aRing is connected
// to.
CxmpInline uint32_t
CxpTrNodeDegree(CxtTr *aTr, CxtTrRing aRing)
{
    uint32_t rVal;
    CxtTrRing ring;

    rVal = 1;
    CxmQriOthersForeach(ring, aTr->trrs, aRing, link)
	{
	    rVal++;
	}

    return rVal;
}

// Calculate the number of edges between two nodes.  A distance of 0 means that
// there is no path between the two nodes.
static uint32_t
CxpTrNodeDistance(CxtTr *aTr, CxtTrRing aRing, CxtTrNode aOther,
		  uint32_t aDistance)
{
    uint32_t rVal;
    CxtTrRing ring;

    if (CxTrRingNodeGet(aTr, aRing) == aOther)
    {
	rVal = aDistance;
	goto RETURN;
    }

    CxmQriOthersForeach(ring, aTr->trrs, aRing, link)
	{
	    if ((rVal = CxpTrNodeDistance(aTr, CxTrRingOtherGet(aTr, ring),
					  aOther, aDistance + 1)) != 0)
	    {
		goto RETURN;
	    }
	}

    rVal = 0;
    RETURN:
    return rVal;
}

CxtTrNode
CxTrNodeNew(CxtTr *aTr)
{
    return CxpTrNodeAlloc(aTr);
}

void
CxTrNodeDelete(CxtTr *aTr, CxtTrNode aNode)
{
    CxmDassert(CxpTrNodeValidate(aTr, aNode));
    CxmAssert(CxmQliFirst(&aTr->trns[aNode].rings) == CxmTrRingNone);

    CxpTrNodeDealloc(aTr, aNode);
}

uint32_t
CxTrNodeTaxonNumGet(CxtTr *aTr, CxtTrNode aNode)
{
    CxmDassert(CxpTrNodeValidate(aTr, aNode));

    return aTr->trns[aNode].taxonNum;
}

void
CxTrNodeTaxonNumSet(CxtTr *aTr, CxtTrNode aNode,
		    uint32_t aTaxonNum)
{
    CxmDassert(CxpTrNodeValidate(aTr, aNode));

    aTr->trns[aNode].taxonNum = aTaxonNum;

    aTr->modified = true;
}

CxtTrRing
CxTrNodeRingGet(CxtTr *aTr, CxtTrNode aNode)
{
    CxtTrRing rVal;
    CxtTrn *trn;

    CxmDassert(CxpTrNodeValidate(aTr, aNode));

    trn = &aTr->trns[aNode];

    rVal = CxmQliFirst(&trn->rings);

    return rVal;
}

void *
CxTrNodeAuxGet(CxtTr *aTr, CxtTrNode aNode)
{
    CxmDassert(CxpTrNodeValidate(aTr, aNode));

    return aTr->trns[aNode].u.aux;
}

void
CxTrNodeAuxSet(CxtTr *aTr, CxtTrNode aNode, void *aAux)
{
    CxmDassert(CxpTrNodeValidate(aTr, aNode));

    aTr->trns[aNode].u.aux = aAux;
}

uint32_t
CxTrNodeDegree(CxtTr *aTr, CxtTrNode aNode)
{
    uint32_t rVal;
    CxtTrRing ring;

    CxmDassert(CxpTrNodeValidate(aTr, aNode));

    ring = CxmQliFirst(&aTr->trns[aNode].rings);
    if (ring != CxmTrRingNone)
    {
	rVal = CxpTrNodeDegree(aTr, ring);
    }
    else
    {
	rVal = 0;
    }

    return rVal;
}

// Set the head of the ring for aNode to aRing.
void
CxTrNodeRingSet(CxtTr *aTr, CxtTrNode aNode, CxtTrRing aRing)
{
    CxmAssert(CxTrRingNodeGet(aTr, aRing) == aNode);

    CxmQliFirst(&aTr->trns[aNode].rings) = aRing;
}

uint32_t
CxTrNodeDistance(CxtTr *aTr, CxtTrNode aNode, CxtTrNode aOther)
{
    uint32_t rVal;
    CxtTrn *trn;
    CxtTrRing ring;

    CxmDassert(CxpTrNodeValidate(aTr, aNode));
    CxmDassert(CxpTrNodeValidate(aTr, aOther));
    CxmAssert(aNode != aOther);

    trn = &aTr->trns[aNode];

    ring = CxmQliFirst(&trn->rings);
    if (ring != CxmTrRingNone)
    {
	CxmQliForeach(ring, &trn->rings, aTr->trrs, link)
	    {
		if ((rVal = CxpTrNodeDistance(aTr, CxTrRingOtherGet(aTr, ring),
					      aOther, 1)) != 0)
		{
		    break;
		}
	    }
    }
    else
    {
	rVal = 0;
    }

    return rVal;
}

//==============================================================================

// tr.

// Initialize everything except trns and sparetrns.
CxmpInline void
CxpTrNew(CxtTr *aTr)
{
    aTr->aux = NULL;
    aTr->modified = false;
    aTr->base = CxmTrNodeNone;
    aTr->ntaxa = 0;
    aTr->nedges = 0;
    aTr->bedges = NULL;
    aTr->nbedgesA = 0;
    aTr->nbedgesB = 0;
    aTr->trt = NULL;
    aTr->trtused = 0;
    aTr->trns = NULL;
    aTr->ntrns = 0;
    aTr->sparetrns = CxmTrNodeNone;
    aTr->tres = NULL;
    aTr->ntres = 0;
    aTr->sparetres = CxmTrEdgeNone;
    aTr->trrs = NULL;
    aTr->held = NULL;
    aTr->heldlen = 0;
    aTr->nheld = 0;

#ifdef CxmDebug
    aTr->magic = CxmTrMagic;
#endif
}

// Recursively traverse the tree, count the number of taxa, and find the lowest
// numbered taxon.
static CxtTrNode
CxpTrLowestRecurse(CxtTr *aTr, CxtTrRing aRing, uint32_t *rNtaxa,
		   uint32_t *rNedges, CxtTrNode aRoot)
{
    CxtTrNode rVal, node, root, troot;
    CxtTrRing ring;
    CxtTrn *trn;

    node = CxTrRingNodeGet(aTr, aRing);
    trn = &aTr->trns[node];

    if (trn->taxonNum != CxmTrNodeTaxonNone)
    {
	// Leaf node.
	(*rNtaxa)++;
    }

    if (trn->taxonNum != CxmTrNodeTaxonNone
	&& (aRoot == CxmTrNodeNone
	    || trn->taxonNum < aTr->trns[aRoot].taxonNum))
    {
	rVal = node;
	root = node;
    }
    else
    {
	rVal = CxmTrNodeNone;
	root = aRoot;
    }

    // Iterate over neighbors.
    CxmQriOthersForeach(ring, aTr->trrs, aRing, link)
	{
	    // Count edge.
	    (*rNedges)++;

	    troot = CxpTrLowestRecurse(aTr, CxTrRingOtherGet(aTr, ring),
				       rNtaxa, rNedges, root);
	    if (troot != CxmTrNodeNone)
	    {
		rVal = troot;
		root = troot;
	    }
	}

    return rVal;
}

// Recursively traverse the tree, count the number of taxa, and find the lowest
// numbered taxon.
static CxtTrNode
CxpTrLowest(CxtTr *aTr, CxtTrNode aNode, uint32_t *rNtaxa,
	    uint32_t *rNedges)
{
    CxtTrNode rVal, root, troot;
    CxtTrRing ring;
    CxtTrn *trn;

    CxmAssert(aNode != CxmTrNodeNone);

    trn = &aTr->trns[aNode];

    if (trn->taxonNum != CxmTrNodeTaxonNone)
    {
	// Leaf node.
	(*rNtaxa)++;
    }

    if (trn->taxonNum != CxmTrNodeTaxonNone)
    {
	rVal = aNode;
	root = aNode;
    }
    else
    {
	rVal = CxmTrNodeNone;
	root = CxmTrNodeNone;
    }

    // Iterate over neighbors.
    CxmQliForeach(ring, &trn->rings, aTr->trrs, link)
	{
	    // Count edge.
	    (*rNedges)++;

	    troot = CxpTrLowestRecurse(aTr, CxTrRingOtherGet(aTr, ring),
				       rNtaxa, rNedges, root);
	    if (troot != CxmTrNodeNone)
	    {
		rVal = troot;
		root = troot;
	    }
	}

    return rVal;
}

#ifdef CxmDebug
// Validate a tree.
bool
CxpTrValidate(CxtTr *aTr)
{
    uint32_t i, ntaxa, nedges;

    CxmCheckPtr(aTr);
    CxmAssert(aTr->magic == CxmTrMagic);
    CxmAssert(aTr->modified == false);

    ntaxa = 0;
    nedges = 0;
    if (aTr->base != CxmTrNodeNone)
    {
	CxpTrLowest(aTr, aTr->base, &ntaxa, &nedges);
    }
    CxmAssert(aTr->ntaxa == ntaxa);
    CxmAssert(aTr->nedges == nedges);

    // Iterate over trns and do some basic sanity checks.
    for (i = 0; i < aTr->ntrns; i++)
    {
	if (aTr->trns[i].magic == CxmTrnMagic)
	{
	    CxpTrNodeValidate(aTr, (CxtTrNode) i);
	}
	else
	{
	    // Make sure there are no valid trn's in the free list.
	    CxmAssert(aTr->trns[i].u.link == CxmTrNodeNone
		      || aTr->trns[aTr->trns[i].u.link].magic
		      != CxmTrnMagic);
	}
    }

    CxmAssert(aTr->sparetrns == CxmTrNodeNone
	      || aTr->trns[aTr->sparetrns].magic != CxmTrnMagic);

    return true;
}
#endif

static void
CxpTrNtaxaNedgesUpdate(CxtTr *aTr)
{
    uint32_t ntaxa, nedges;

    // Update ntaxa and nedges.
    ntaxa = 0;
    nedges = 0;
    if (aTr->base != CxmTrNodeNone)
    {
	CxpTrLowest(aTr, aTr->base, &ntaxa, &nedges);
    }

    aTr->ntaxa = ntaxa;
    aTr->nedges = nedges;
}

static void
CxpTrBisectionEdgeListGenRecurse(CxtTr *aTr, CxtTrRing aRing,
				 CxtTrEdge *arEdges,
				 uint32_t *arNedges)
{
    CxtTrRing ring;

    CxmQriOthersForeach(ring, aTr->trrs, aRing, link)
	{
	    // Add edge to list.
	    arEdges[*arNedges] = CxTrRingEdgeGet(aTr, ring);
	    (*arNedges)++;

	    // Recurse into neighbor subtree.
	    CxpTrBisectionEdgeListGenRecurse(aTr,
					     CxTrRingOtherGet(aTr, ring),
					     arEdges, arNedges);
	}
}

// Pretend that the tree is bisected at the edge that contains aRing.
// Construct a list of edges that are in the subtree that contains aRing.
//
// The first element in the list is always the edge that is adjacent to the
// bisection.  This facilitates recognition of reconnections that would reverse
// bisection.
//
// If the list is empty (bisection adjacent to a leaf node), return the single
// node, so that it can be accessed directly (there's no edge logically attached
// to it).
CxmpInline CxtTrNode
CxpTrBisectionEdgeListGen(CxtTr *aTr, CxtTrRing aRing,
			  CxtTrEdge *arEdges, uint32_t *arNedges)
{
    CxtTrNode rVal;
    CxtTrRing ring;

    // Initialize the length of the list before recursing.
    *arNedges = 0;

    switch (CxpTrNodeDegree(aTr, aRing))
    {
	case 1:
	{
	    // A subtree that is composed of a single node has no edges.  Add a
	    // single entry to the list, and return the node.
	    arEdges[0] = CxmTrEdgeNone;
	    (*arNedges)++;
	    rVal = CxTrRingNodeGet(aTr, aRing);
	    break;
	}
	case 2:
	{
	    // A tree should never have nodes of degree 2.
	    CxmNotReached();
	}
	case 3:
	{
	    // Take care to add only one of the edges that is connected to the
	    // node, since from the perspective of TBR, the node does not exist.
	    // (A node of degree 2 is a superfluous node.)

	    // First edge.
	    ring = CxmQriNext(aTr->trrs, aRing, link);
	    arEdges[0] = CxTrRingEdgeGet(aTr, ring);
	    (*arNedges)++;

	    // First subtree.
	    CxpTrBisectionEdgeListGenRecurse(aTr,
					     CxTrRingOtherGet(aTr,
							       ring),
					     arEdges, arNedges);

	    // Second subtree.
	    ring = CxmQriNext(aTr->trrs, ring, link);
	    CxpTrBisectionEdgeListGenRecurse(aTr,
					     CxTrRingOtherGet(aTr,
							       ring),
					     arEdges, arNedges);

	    rVal = CxmTrNodeNone;
	    break;
	}
	default:
	{
	    // Add all edges in the subtree.  Removing the bisection edge still
	    // leaves enough edges attached to the node for the node to have
	    // relevance.
	    CxmQriOthersForeach(ring, aTr->trrs, aRing, link)
		{
		    // Add edge to list.
		    arEdges[*arNedges] = CxTrRingEdgeGet(aTr, ring);
		    (*arNedges)++;

		    CxpTrBisectionEdgeListGenRecurse(aTr,
						     CxTrRingOtherGet(aTr,
								      ring),
						     arEdges, arNedges);
		}

	    rVal = CxmTrNodeNone;
	    break;
	}
    }

    return rVal;
}

// Generate lists of edges in each half of a logical bisection at edge
// aBisect.
CxmpInline void
CxpTrBedgesGen(CxtTr *aTr, CxtTrEdge aBisect, CxtTrNode *rNodeA,
	       CxtTrNode *rNodeB)
{
    CxtTrNode nodeA, nodeB;

    CxmDassert(CxpTrEdgeValidate(aTr, aBisect));

    nodeA = CxpTrBisectionEdgeListGen(aTr,
				      CxTrEdgeRingGet(aTr, aBisect, 0),
				      aTr->bedges, &aTr->nbedgesA);
    nodeB = CxpTrBisectionEdgeListGen(aTr,
				      CxTrEdgeRingGet(aTr, aBisect, 1),
				      &aTr->bedges[aTr->nbedgesA],
				      &aTr->nbedgesB);

    if (rNodeA != NULL)
    {
	*rNodeA = nodeA;
    }
    if (rNodeB != NULL)
    {
	*rNodeB = nodeB;
    }
}

static void
CxpTrtBisectEdgeUpdateRecurse(CxtTr *aTr, CxtTrRing aRing,
			      uint32_t *arEdgeCount)
{
    CxtTrRing ring;

    CxmQriOthersForeach(ring, aTr->trrs, aRing, link)
	{
	    // Record edge.
	    aTr->trt[*arEdgeCount].bisectEdge = CxTrRingEdgeGet(aTr, ring);
	    (*arEdgeCount)++;

	    // Recurse into neighbor subtree.
	    CxpTrtBisectEdgeUpdateRecurse(aTr,
					  CxTrRingOtherGet(aTr, ring),
					  arEdgeCount);
	}
}

static void
CxpTrtUpdate(CxtTr *aTr, uint32_t aNedgesPrev)
{
    uint32_t i, j, n, offset, degreeA, degreeB;

    CxmAssert(aTr->modified == false);

    // Allocate/reallocate/deallocate trt.
    if (aTr->trt == NULL)
    {
	// Allocate trt.
	aTr->trt = (CxtTrt *) CxmMalloc(sizeof(CxtTrt)
					* (aTr->nedges + 1));
    }
    else if (aTr->nedges != aNedgesPrev)
    {
	// Reallocate trt.  There is never a need to deallocate trt here,
	// since trt contains one extra element.
	aTr->trt = (CxtTrt *) CxmRealloc(aTr->trt,
					 sizeof(CxtTrt)
					 * (aTr->nedges + 1));
    }

    // Recursively traverse the tree, and initialize trt->bisectEdge along the
    // way.
    if (aTr->nedges > 0)
    {
	uint32_t edgeCount;
	CxtTrn *trn;
	CxtTrRing ring;

	CxmAssert(aTr->base != CxmTrNodeNone);

	edgeCount = 0;
	trn = &aTr->trns[aTr->base];
	CxmQliForeach(ring, &trn->rings, aTr->trrs, link)
	    {
		// Record edge.
		aTr->trt[edgeCount].bisectEdge = CxTrRingEdgeGet(aTr, ring);
		edgeCount++;

		CxpTrtBisectEdgeUpdateRecurse(aTr,
					      CxTrRingOtherGet(aTr, ring),
					      &edgeCount);
	    }
	CxmAssert(edgeCount == aTr->nedges);
    }

    // Iteratively fill in trt.
    for (i = j = offset = 0; i < aTr->nedges; i++)
    {
	// Record offset.
	aTr->trt[j].offset = offset;

	// Update offset.
	CxpTrBedgesGen(aTr, aTr->trt[i].bisectEdge, NULL, NULL);
	n = (aTr->nbedgesA * aTr->nbedgesB);
	// TBR can only be reversed for an edge that is attached to two nodes of
	// degree 3.
	degreeA = CxpTrNodeDegree(aTr, CxTrEdgeRingGet(aTr,
						       aTr->trt[i].bisectEdge,
						       0));
	degreeB = CxpTrNodeDegree(aTr, CxTrEdgeRingGet(aTr,
						       aTr->trt[i].bisectEdge,
						       1));
	if ((degreeA == 1 || degreeA == 3) && (degreeB == 1 || degreeB == 3))
	{
	    n--;
	}
	if (n != 0)
	{
	    aTr->trt[j].bisectEdge = aTr->trt[i].bisectEdge;
	    offset += n;
	    j++;
	}
    }
    aTr->trt[j].offset = offset;

    // It may be that not all bisections result in neighbors, so the table may
    // not be full.  Keep track of the number of valid elements (not counting
    // the trailing one that stores the total number of TBR neighbors).
    aTr->trtused = j;
}

static void
CxpTrBedgesUpdate(CxtTr *aTr, uint32_t aNedgesPrev)
{
    // Allocate/reallocate/deallocate bedges.  To keep things simple, allocate
    // as big an array as there are edges, even though not quite that many are
    // ever used.
    if (aTr->bedges == NULL)
    {
	// Allocate bedges.
	aTr->bedges
	    = (CxtTrEdge *) CxmMalloc(sizeof(CxtTrEdge) * aTr->nedges);
    }
    else if (aTr->nedges != aNedgesPrev)
    {
	if (aTr->nedges > 0)
	{
	    // Reallocate bedges.
	    aTr->bedges = (CxtTrEdge *) CxmRealloc(aTr->bedges,
						   sizeof(CxtTrEdge)
						   * aTr->nedges);
	}
	else
	{
	    // Deallocate bedges.
	    CxmFree(aTr->bedges);
	    aTr->bedges = NULL;
	}
    }

    // Clear nbedges_[ab].
    aTr->nbedgesA = 0;
    aTr->nbedgesB = 0;
}

CxmpInline void
CxpTrUpdate(CxtTr *aTr)
{
    if (aTr->modified)
    {
	uint32_t nedgesPrev;

	// Store nedges before updating.
	nedgesPrev = aTr->nedges;

	// Update ntaxa and nedges.
	CxpTrNtaxaNedgesUpdate(aTr);

	// Reset the modified flag.
	aTr->modified = false;

	// Update bedges and trt.
	CxpTrBedgesUpdate(aTr, nedgesPrev);
	CxpTrtUpdate(aTr, nedgesPrev);

	// Clear held trees.
	aTr->nheld = 0;
    }
}

static void
CxpTrMpRingFinish(CxtTr *aTr, CxtTrRing aRing)
{
    CxtTrr *trr;

    trr = &aTr->trrs[aRing];

    if (trr->ps != NULL)
    {
	CxpTrPsDelete(aTr, trr->ps);
	trr->ps = NULL;
    }
}

static void
CxpTrMpFinishRecurse(CxtTr *aTr, CxtTrRing aRing)
{
    CxtTrRing ring;
    CxtTre *tre;

    // Clean up aRing.
    CxpTrMpRingFinish(aTr, aRing);

    // Recurse into subtrees.
    CxmQriOthersForeach(ring, aTr->trrs, aRing, link)
	{
	    // Clean up edge before recursing.
	    tre = &aTr->tres[CxTrRingEdgeGet(aTr, ring)];
	    if (tre->ps != NULL)
	    {
		CxpTrPsDelete(aTr, tre->ps);
		tre->ps = NULL;
	    }

	    // Clean up ring.
	    CxpTrMpRingFinish(aTr, ring);

	    // Recurse.
	    CxpTrMpFinishRecurse(aTr, CxTrRingOtherGet(aTr, ring));
	}
}

#ifdef CxmCpuIa32
CxmpInline void
CxpTrMpIa32Pscore(CxtTr *aTr, CxtTrPs *aP, CxtTrPs *aA,
		  CxtTrPs *aB)
{
    uint32_t curlimit, i, nbytes, ns;
    CxtTrc *charsP, *charsA, *charsB;

    // Calculate node score.
    ns = 0;

    // Calculate partial Fitch parsimony scores for each character.
    charsP = aP->chars;
    charsA = aA->chars;
    charsB = aB->chars;

    nbytes = (aP->nchars >> 1);

    // Initialize SSE2 registers.
    {
	static const unsigned char low[] __attribute__ ((aligned (16))) =
	    "\x0f\x0f\x0f\x0f\x0f\x0f\x0f\x0f"
	    "\x0f\x0f\x0f\x0f\x0f\x0f\x0f\x0f";

	asm volatile (
	    // Clear pns.
	    "pxor %%xmm4, %%xmm4;"

	    // Fill xmm5 with masks for the least significant four bits of each
	    // byte.
	    "movdqa %[low], %%xmm5;"

	    // Fill xmm6 with masks for the most significant four bits of each
	    // byte.
	    "pcmpeqb %%xmm6, %%xmm6;"
	    "pxor %%xmm5, %%xmm6;"

	    // Fill xmm7 with 16 1's.
	    "pxor %%xmm7, %%xmm7;"
	    "pcmpeqb %%xmm0, %%xmm0;"
	    "psubb %%xmm0, %%xmm7;"
	    :
	    : [low] "m" (*low)
	    : "%xmm4", "%xmm5", "%xmm6", "%xmm7"
	    );
    }

    // The inner loop can be run a maximum of 127 times before the partial node
    // score results (stored in %xmm4) are added to ns (otherwise, overflow
    // could occur).  Therefore, the outer loop calculates the upper bound for
    // the inner loop, thereby avoiding extra computation in the inner loop.
    curlimit = 127 * 16;
    if (curlimit > nbytes)
    {
	curlimit = nbytes;
    }
    for (i = 0;;)
    {
	// Use SSE2 to evaluate the characters.  This loop handles 32 characters
	// per iteration.
	for (; i < curlimit; i += 16)
	{
	    asm volatile (
		// Read character data, and'ing and or'ing them together.
		//
		// a = *charsA;
		// b = *charsB;
		// p = a & b;
		// d = a | b;
		"movdqa %[a], %%xmm0;"
		"movdqa %%xmm0, %%xmm1;"
		"por %[b], %%xmm0;" // xmm0 contains d.
		"pand %[b], %%xmm1;" // xmm1 contains p.

		//==============================================================
		// Most significant bits.

		// Create bitmasks according to whether the character state sets
		// are empty.
		//
		// c = p ? 0x00 : 0xff;
		// e = (c & d);
		// s = c ? 0 : 1;
		// p = (p | e);
		// ns += s;
		"pxor %%xmm2, %%xmm2;"
		"movdqa %%xmm1, %%xmm3;"
		"pand %%xmm6, %%xmm3;" // Mask out unwanted bits of p.
		"pcmpeqb %%xmm3, %%xmm2;" // xmm2 contains c.
		"movdqa %%xmm0, %%xmm3;"
		"pand %%xmm6, %%xmm3;" // Mask out unwanted bits of d.
		"pand %%xmm2, %%xmm3;" // xmm3 contains e.
		"por %%xmm3, %%xmm1;" // Update p.
		"pand %%xmm7, %%xmm2;" // xmm2 contains s.
		"paddusb %%xmm2, %%xmm4;" // Update ns (add s).

		//==============================================================
		// Least significant bits.

		// Create bitmasks according to whether the character state sets
		// are empty.
		//
		// c = p ? 0x00 : 0xff;
		// e = (c & d);
		// s = c ? 0 : 1;
		// p = (p | e);
		// ns += s;
		"pxor %%xmm2, %%xmm2;"
		"movdqa %%xmm1, %%xmm3;"
		"pand %%xmm5, %%xmm3;" // Mask out unwanted bits of p.
		"pcmpeqb %%xmm3, %%xmm2;" // xmm2 contains c.
		"pand %%xmm5, %%xmm0;" // Mask out unwanted bits of d.
		"pand %%xmm2, %%xmm0;" // xmm0 contains e.
		"por %%xmm0, %%xmm1;" // Update p.
		"pand %%xmm7, %%xmm2;" // xmm2 contains s.
		"paddusb %%xmm2, %%xmm4;" // Update ns (add s).

		// Store results.
		//
		// *charsP = p;
		"movdqa %%xmm1, %[p];"
		: [p] "=m" (charsP[i])
		: [a] "m" (charsA[i]), [b] "m" (charsB[i])
		: "memory"
		);
	}

	// Update ns and reset pns.
	{
	    uint32_t pns;

	    asm volatile (
		// Sum the upper 8 bytes and lower 8 bytes separately (that's
		// what psadbw does).
		"pxor %%xmm0, %%xmm0;"
		"psadbw %%xmm0, %%xmm4;"

		// Combine the results of psadbw.
		"movdqa %%xmm4, %%xmm3;"
		"punpckhqdq %%xmm0, %%xmm4;"
		"paddq %%xmm3, %%xmm4;"

		// Store the result.
		"movd %%xmm4, %[pns];"
		"pxor %%xmm4, %%xmm4;"
		: [pns] "=r" (pns)
		:
		: "memory"
		);

	    ns += pns;
	}

	// Break out of the loop if the bound for the inner loop was the maximum
	// possible.
	if (curlimit == nbytes)
	{
	    break;
	}
	// Update the bound for the inner loop, taking care not to exceed the
	// maximum possible bound.
	curlimit += 127 * 16;
	if (curlimit > nbytes)
	{
	    curlimit = nbytes;
	}
    }

    aP->nodeScore = ns;
}
#endif

static void
CxpTrMpCPscore(CxtTr *aTr, CxtTrPs *aP, CxtTrPs *aA,
	       CxtTrPs *aB)
{
    uint32_t i, nwords, ns, a, b, m, r, un;
    uint32_t *charsP, *charsA, *charsB;
    static const uint32_t bitsTable[] =
    {
	2, 1, -1, -1, -1, -1, -1, -1,
	-1, -1, -1, -1, -1, -1, -1, -1,
	1, 0
    };

    // Calculate node score.
    ns = 0;

#define CxmTrMpCPscoreInner()						\
    a = charsA[i];							\
    b = charsB[i];							\
									\
    /* Get 1's in the least significant bits of state sets that are	\
     * non-empty after the intersection operation. */			\
    r = m = a & b;							\
    m |= (m >> 1);							\
    m |= (m >> 1);							\
    m |= (m >> 1);							\
									\
    /* Mask out garbage. */						\
    m &= 0x11111111;							\
									\
    /* Count up changes. */						\
    ns += bitsTable[m & 0xff]						\
	+ bitsTable[(m >> 8) & 0xff]					\
	+ bitsTable[(m >> 16) & 0xff]					\
	+ bitsTable[(m >> 24) & 0xff];					\
									\
    /* Propagate 1's to make a bit mask. */				\
    m |= (m << 1);							\
    m |= (m << 1);							\
    m |= (m << 1);							\
									\
    /* Combine results of intersection and union operations. */		\
    r &= m;								\
    un = a | b;								\
    m = (~m);								\
    un &= m;								\
    r |= un;								\
									\
    /* Store result. */							\
    charsP[i] = r;							\
									\
    i++;

    // Calculate preliminary Fitch parsimony scores for each character.
    charsP = (uint32_t *) aP->chars;
    charsA = (uint32_t *) aA->chars;
    charsB = (uint32_t *) aB->chars;
    for (i = 0, nwords = (aP->nchars >> 3); i < nwords;)
    {
	CxmTrMpCPscoreInner();
	CxmTrMpCPscoreInner();
	CxmTrMpCPscoreInner();
	CxmTrMpCPscoreInner();
    }
#undef CxmTrMpCPscoreInner

    aP->nodeScore = ns;
}

// Unconditionally calculate the partial score for aP, using aA and aB as
// children.
CxmpInline void
CxpTrMpPscore(CxtTr *aTr, CxtTrPs *aP, CxtTrPs *aA, CxtTrPs *aB)
{
    // Reset this node's parent pointer, to keep the parent from using an
    // invalid cached value.
    aP->parent = NULL;

    // Calculate sum of subtree scores.
    aP->subtreesScore
	= aA->subtreesScore + aA->nodeScore
	+ aB->subtreesScore + aB->nodeScore;

#ifdef CxmCpuIa32
    if (CxgIa32UseSse2)
    {
	CxpTrMpIa32Pscore(aTr, aP, aA, aB);
    }
    else
#endif
    {
	CxpTrMpCPscore(aTr, aP, aA, aB);
    }
}

// The sole purpose of this function is to assure that the contents of
// CxpTrMpPscore() are not inlined in CxpTrMpCachePscore().  Most of the time,
// the cache should be usable, so the actual scoring code doesn't usually get
// called.
static void
CxpTrNoInlineMpPscore(CxtTr *aTr, CxtTrPs *aP, CxtTrPs *aA,
		      CxtTrPs *aB)
{
    CxpTrMpPscore(aTr, aP, aA, aB);
}

// Calculate the partial score for aP, using aA and aB as children.  However, do
// some extra bookkeeping in order to be able to cache the results, and later
// recognize that precisely the same calculation was cached.
CxmpInline void
CxpTrMpCachePscore(CxtTr *aTr, CxtTrPs *aP, CxtTrPs *aA,
		   CxtTrPs *aB)
{
//#define CxmTrMpCachePscoreValidate
#ifdef CxmTrMpCachePscoreValidate
    bool cached;
    uint32_t cachedNodeScore;
#endif

    CxmCheckPtr(aP);
    CxmCheckPtr(aA);
    CxmCheckPtr(aB);

    // Only calculate the parent's node score if the cached value is invalid.
    if (aA->parent != aP || aB->parent != aP)
#ifdef CxmTrMpCachePscoreValidate
    {
	cached = false;
    }
    else
    {
	cached = true;
	cachedNodeScore = aP->nodeScore;

	if (aP->subtreesScore
	    != (aA->subtreesScore + aA->nodeScore
		+ aB->subtreesScore + aB->nodeScore))
	{
	    fprintf(stderr,
		    "%s:%d:%s(): subtreesScore %u (should be %u)\n",
		    __FILE__, __LINE__, __FUNCTION__,
		    aP->subtreesScore,
		    aA->subtreesScore + aA->nodeScore
		    + aB->subtreesScore + aB->nodeScore);
	    abort();
	}
    }
#endif
    {
	// Set parent pointers, so that cached values may be used in future
	// runs.
	aA->parent = aP;
	aB->parent = aP;

	// Calculate the partial score.
	CxpTrNoInlineMpPscore(aTr, aP, aA, aB);
    }

#ifdef CxmTrMpCachePscoreValidate
    if (cached)
    {
	if (cachedNodeScore != aP->nodeScore)
	{
	    fprintf(stderr, "%s:%d:%s(): nodeScore %u (should be %u)\n",
		    __FILE__, __LINE__, __FUNCTION__,
		    cachedNodeScore, aP->nodeScore);
	    abort();
	}
    }
#endif
}

CxmpInline void
CxpTrMpCacheInvalidate(CxtTr *aTr, CxtTrPs *aPs)
{
    CxmCheckPtr(aPs);

    // Reset this node's parent pointer, to keep the old parent from using an
    // invalid cached value.
    aPs->parent = NULL;
}

#ifdef CxmCpuIa32
CxmpInline uint32_t
CxpTrMpIa32Fscore(CxtTr *aTr, CxtTrPs *aA, CxtTrPs *aB,
		  uint32_t aMaxscore)
{
    uint32_t rVal, i, nbytes, pns;
    CxtTrc *charsA, *charsB;

    // Calculate sum of subtree scores.
    rVal
	= aA->subtreesScore + aA->nodeScore
	+ aB->subtreesScore + aB->nodeScore;

    // Calculate partial Fitch parsimony scores for each character.
    charsA = aA->chars;
    charsB = aB->chars;

    // Initialize SSE2 registers.
    {
	static const unsigned char low[] __attribute__ ((aligned (16))) =
	    "\x0f\x0f\x0f\x0f\x0f\x0f\x0f\x0f"
	    "\x0f\x0f\x0f\x0f\x0f\x0f\x0f\x0f";

	asm volatile (
	    // Clear pns.
	    "pxor %%xmm4, %%xmm4;"

	    // Fill xmm5 with masks for the least significant four bits of each
	    // byte.
	    "movdqa %[low], %%xmm5;"

	    // Fill xmm6 with masks for the most significant four bits of each
	    // byte.
	    "pcmpeqb %%xmm6, %%xmm6;"
	    "pxor %%xmm5, %%xmm6;"

	    // Fill xmm7 with 16 1's.
	    "pxor %%xmm7, %%xmm7;"
	    "pcmpeqb %%xmm0, %%xmm0;"
	    "psubb %%xmm0, %%xmm7;"
	    :
	    : [low] "m" (*low)
	    : "%xmm4", "%xmm5", "%xmm6", "%xmm7"
	    );
    }

    // Use SSE2 to evaluate the characters.  This loop handles 32 characters per
    // iteration.
    for (i = 0, nbytes = (aA->nchars >> 1); i < nbytes; i += 16)
    {
	asm volatile (
	    // Read character data, and'ing and or'ing them together.
	    //
	    // a = *charsA;
	    // b = *charsB;
	    // p = a & b;
	    "movdqa %[a], %%xmm1;"
	    "pand %[b], %%xmm1;" // xmm1 contains p.

	    //==================================================================
	    // Most significant bits.

	    // Create bitmasks according to whether the character state sets are
	    // empty.
	    //
	    // c = p ? 0x00 : 0xff;
	    // s = c ? 0 : 1;
	    // rVal += s;
	    "pxor %%xmm2, %%xmm2;"
	    "movdqa %%xmm1, %%xmm3;"
	    "pand %%xmm6, %%xmm3;" // Mask out unwanted bits of p.
	    "pcmpeqb %%xmm3, %%xmm2;" // xmm2 contains c.
	    "pand %%xmm7, %%xmm2;" // xmm2 contains s.
	    "paddusb %%xmm2, %%xmm4;" // Update rVal (add s).

	    //==================================================================
	    // Least significant bits.

	    // Create bitmasks according to whether the character state sets are
	    // empty.
	    //
	    // c = p ? 0x00 : 0xff;
	    // s = c ? 0 : 1;
	    // rVal += s;
	    "pxor %%xmm2, %%xmm2;"
	    "pand %%xmm5, %%xmm1;" // Mask out unwanted bits of p.
	    "pcmpeqb %%xmm1, %%xmm2;" // xmm2 contains c.
	    "pand %%xmm7, %%xmm2;" // xmm2 contains s.
	    "paddusb %%xmm2, %%xmm4;" // Update rVal (add s).

	    // Sum the upper 8 bytes and lower 8 bytes separately (there's no
	    // choice in the matter -- that's what psadbw does).
	    "pxor %%xmm0, %%xmm0;"
	    "psadbw %%xmm0, %%xmm4;"

	    // Combine the results of psadbw.
	    "movdqa %%xmm4, %%xmm3;"
	    "punpckhqdq %%xmm0, %%xmm4;"
	    "paddq %%xmm3, %%xmm4;"

	    // Store the result.
	    "movd %%xmm4, %[pns];"
	    "pxor %%xmm4, %%xmm4;"

	    : [pns] "=r" (pns)
	    : [a] "m" (charsA[i]), [b] "m" (charsB[i])
	    : "memory"
	    );

	// Update rVal and terminate if the max score was exceeded.
	rVal += pns;
	if (rVal > aMaxscore)
	{
	    rVal = UINT_MAX;
	    break;
	}
    }

    return rVal;
}
#endif

static uint32_t
CxpTrMpCFscore(CxtTr *aTr, CxtTrPs *aA, CxtTrPs *aB,
	       uint32_t aMaxscore)
{
    uint32_t rVal, i, nwords, a, b, m;
    uint32_t *charsA, *charsB;
    static const uint32_t bitsTable[] =
    {
	2, 1, -1, -1, -1, -1, -1, -1,
	-1, -1, -1, -1, -1, -1, -1, -1,
	1, 0
    };

    // Calculate sum of subtree scores.
    rVal
	= aA->subtreesScore + aA->nodeScore
	+ aB->subtreesScore + aB->nodeScore;

#define CxmTrMpCFscoreInner()						\
    a = charsA[i];							\
    b = charsB[i];							\
									\
    /* Get 1's in the least significant bits of state sets that are	\
     * non-empty after the intersection operation. */			\
    m = a & b;								\
    m |= (m >> 1);							\
    m |= (m >> 1);							\
    m |= (m >> 1);							\
									\
    /* Mask out garbage. */						\
    m &= 0x11111111;							\
									\
    /* Count up changes. */						\
    rVal += bitsTable[m & 0xff]						\
	+ bitsTable[(m >> 8) & 0xff]					\
	+ bitsTable[(m >> 16) & 0xff]					\
	+ bitsTable[(m >> 24) & 0xff];					\
									\
    if (rVal > aMaxscore)						\
    {									\
	rVal = UINT_MAX;						\
	break;								\
    }									\
									\
    i++;

    // Calculate partial Fitch parsimony scores for each character.
    charsA = (uint32_t *) aA->chars;
    charsB = (uint32_t *) aB->chars;
    for (i = 0, nwords = (aA->nchars >> 3); i < nwords;)
    {
	CxmTrMpCFscoreInner();
	CxmTrMpCFscoreInner();
	CxmTrMpCFscoreInner();
	CxmTrMpCFscoreInner();
    }
#undef CxmTrMpCFscoreInner

    return rVal;
}

// Unconditionally calculate the final score of a tree, using aA and aB as
// children.
CxmpInline uint32_t
CxpTrMpFscore(CxtTr *aTr, CxtTrPs *aA, CxtTrPs *aB,
	      uint32_t aMaxscore)
{
    uint32_t rVal;

#ifdef CxmCpuIa32
    if (CxgIa32UseSse2)
    {
	rVal = CxpTrMpIa32Fscore(aTr, aA, aB, aMaxscore);
    }
    else
#endif
    {
	rVal = CxpTrMpCFscore(aTr, aA, aB, aMaxscore);
    }

    return rVal;
}

static CxtTrPs *
CxpTrMpScoreRecurse(CxtTr *aTr, CxtTrRing aRing, CxtTrEdge aBisect)
{
    CxtTrPs *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    uint32_t degree;
    bool adjacent;
    CxtTrRing ring;

    // Get the degree of the node.  Don't count the bisection edge (only an
    // issue if this node is adjacent to the bisection).
    degree = 1;
    adjacent = false;
    CxmQriOthersForeach(ring, aTr->trrs, aRing, link)
	{
	    if (CxTrRingEdgeGet(aTr, ring) != aBisect)
	    {
		degree++;
	    }
	    else
	    {
		adjacent = true;
	    }
	}

    switch (degree)
    {
	case 1:
	{
	    // Leaf node.  Do nothing.
	    rVal = aTr->trrs[aRing].ps;
	    break;
	}
	case 2:
	{
	    // This is a trifurcating node that is adjacent to the bisection.
	    // Return the child node's ps, since this node's ps is
	    // irrelevant.
	    CxmAssert(adjacent);

	    // Clear the cache for the view that is being bypassed.  This is
	    // critical to correctness of the caching machinery, since each view
	    // should never be claimed as the parent of more than two other
	    // views.
	    CxpTrMpCacheInvalidate(aTr, aTr->trrs[aRing].ps);

	    // Get the ring element that connects to the other portion of the
	    // subtree on this side of the bisection.
	    CxmQriOthersForeach(ring, aTr->trrs, aRing, link)
		{
		    if (CxTrRingEdgeGet(aTr, ring) != aBisect)
		    {
			rVal
			    = CxpTrMpScoreRecurse(aTr,
						  CxTrRingOtherGet(aTr, ring),
						  aBisect);
			break;
		    }
		}
	    break;
	}
	case 3:
	{
	    if (adjacent == false)
	    {
		CxtTrPs *psA, *psB;

		// This is a normal trifurcating node.  This is the common case,
		// and is handled separately from the code below for performance
		// reasons.

		// Recursively calculate partial scores for the subtrees.
		ring = CxmQriNext(aTr->trrs, aRing, link);
		psA = CxpTrMpScoreRecurse(aTr,
					  CxTrRingOtherGet(aTr, ring),
					  aBisect);

		ring = CxmQriNext(aTr->trrs, ring, link);
		psB = CxpTrMpScoreRecurse(aTr,
					  CxTrRingOtherGet(aTr, ring),
					  aBisect);

		// Calculate the partial score for this node.
		rVal = aTr->trrs[aRing].ps;
		CxpTrMpCachePscore(aTr, rVal, psA, psB);

		break;
	    }
	    // Fall through if this node is adjacent to the bisection.
	}
	default:
	{
	    // This is a multifurcating node.
	    CxmError("XXX Not implemented");
	}
    }

    return rVal;
}

static void
CxpTrMpViewsRecurse(CxtTr *aTr, CxtTrRing aRing, CxtTrPs *aPs,
		    CxtTrEdge aBisect)
{
    uint32_t degree;
    bool adjacent;
    CxtTrRing ring;

    // Get the degree of the node.  Don't count the bisection edge (only an
    // issue if this node is adjacent to the bisection).
    degree = 1;
    adjacent = false;
    CxmQriOthersForeach(ring, aTr->trrs, aRing, link)
	{
	    if (CxTrRingEdgeGet(aTr, ring) != aBisect)
	    {
		degree++;
	    }
	    else
	    {
		adjacent = true;
	    }
	}

    switch (degree)
    {
	case 1:
	{
	    // Leaf node.  Do nothing.
	    break;
	}
	case 2:
	{
	    // This is a trifurcating node that is adjacent to the bisection.
	    // Pass the parent's ps when recursing, since this node's ps is
	    // irrelevant.
	    CxmAssert(adjacent);

	    // Get the ring element that connects to the other portion of the
	    // subtree on this side of the bisection.
	    CxmQriOthersForeach(ring, aTr->trrs, aRing, link)
		{
		    if (CxTrRingEdgeGet(aTr, ring) != aBisect)
		    {
			// Clear the cache for the view that is being bypassed.
			// This is critical to correctness of the caching
			// machinery, since each view should never be claimed as
			// the parent of more than two other views.
			CxpTrMpCacheInvalidate(aTr, aTr->trrs[ring].ps);

			// Recurse.
			CxpTrMpViewsRecurse(aTr,
					    CxTrRingOtherGet(aTr, ring),
					    aPs,
					    aBisect);
			break;
		    }
		}
	    break;
	}
	case 3:
	{
	    if (adjacent == false)
	    {
		CxtTrRing ringA, ringB;
		CxtTrRing ringAOther, ringBOther;
		CxtTrPs *psA, *psB;
		CxtTrPs *psAOther, *psBOther;

		// This is a normal trifurcating node.  This is the common case,
		// and is handled separately from the code below for performance
		// reasons.

		// Get all variables that are necessary for view calculation and
		// recursion.
		ringA = CxmQriNext(aTr->trrs, aRing, link);
		psA = aTr->trrs[ringA].ps;
		ringAOther = CxTrRingOtherGet(aTr, ringA);
		psAOther = aTr->trrs[ringAOther].ps;

		ringB = CxmQriNext(aTr->trrs, ringA, link);
		psB = aTr->trrs[ringB].ps;
		ringBOther = CxTrRingOtherGet(aTr, ringB);
		psBOther = aTr->trrs[ringBOther].ps;

		// Calculate views and edges, and recurse.
		CxpTrMpPscore(aTr, psA, aPs, psBOther);
		CxpTrMpPscore(aTr,
			      aTr->tres[CxTrRingEdgeGet(aTr, ringA)].ps,
			      psA,
			      aTr->trrs[ringAOther].ps);
		CxpTrMpViewsRecurse(aTr, ringAOther, psA, aBisect);

		CxpTrMpPscore(aTr, psB, aPs, psAOther);
		CxpTrMpPscore(aTr,
			      aTr->tres[CxTrRingEdgeGet(aTr, ringB)].ps,
			      psB,
			      aTr->trrs[ringBOther].ps);
		CxpTrMpViewsRecurse(aTr, ringBOther, psB, aBisect);

		break;
	    }
	    // Fall through if this node is adjacent to the bisection.
	}
	default:
	{
	    // This is a multifurcating node.
	    CxmError("XXX Not implemented");
	}
    }
}

// Calculate the partial score for each edge in aEdges.  aEdges[0] must either
// be CxmTrEdgeNone, or the edge connected to the node that is in turn connected
// to the bisection edge.
CxmpInline bool
CxpTrBisectionEdgeListMp(CxtTr *aTr, CxtTrEdge *aEdges,
			 uint32_t aNedges, CxtTrEdge aBisect,
			 uint32_t aMaxscore)
{
    bool rVal;

    if (aEdges[0] != CxmTrEdgeNone)
    {
	CxtTrRing ringA, ringB;
	CxtTrPs *ps, *psA, *psB;

	ringA = CxTrEdgeRingGet(aTr, aEdges[0], 0);
	ringB = CxTrEdgeRingGet(aTr, aEdges[0], 1);

	// Recursively (post-order traversal) calculate the partial score at
	// each node, as viewed from the first edge in aEdges.  This leaves one
	// valid view at each node, which then makes it possible to calculate
	// the rest of the views during a pre-order traversal of the tree.
	psA = CxpTrMpScoreRecurse(aTr, ringA, aBisect);
	psB = CxpTrMpScoreRecurse(aTr, ringB, aBisect);

	// The first edge must be calculated using psA and psB as children,
	// rather than using the ps's at the ends of the edge.  This is because
	// one of the connected nodes is in turn connected to the bisection
	// edge, which means that the node does not have a useful ps.  The first
	// edge is the only one for which this is an issue, so it is handled
	// here.
	ps = aTr->tres[aEdges[0]].ps;
	CxpTrMpPscore(aTr, ps, psA, psB);
	if (ps->subtreesScore + ps->nodeScore > aMaxscore)
	{
	    // Don't bother calculating other views or edge states, since this
	    // subtree exceeds the maximum score that's of interest.
	    rVal = true;
	    goto RETURN;
	}

	// Perform the pre-order traversal, calculating the remaining views that
	// were not calculated by the above post-order traversal, as well as
	// calculating the state sets for the edges along the way.  Take care to
	// pass the appropriate ps's.
	CxpTrMpViewsRecurse(aTr, ringA, psB, aBisect);
	CxpTrMpViewsRecurse(aTr, ringB, psA, aBisect);

#ifdef CxmDebug
	// Validate per-edge partial scores.
	{
	    uint32_t i;

	    for (i = 1; i < aNedges; i++)
	    {
		// All edge partial scores should have the same value, since the
		// location of the root is irrelevant to the score.
		if (aTr->tres[aEdges[i]].ps->subtreesScore
		    + aTr->tres[aEdges[i]].ps->nodeScore
		    != aTr->tres[aEdges[0]].ps->subtreesScore
		    + aTr->tres[aEdges[0]].ps->nodeScore)
		{
		    fprintf(stderr,
			    "%s:%d:%s(): Expected %u (%u + %u),"
			    " got %u (%u + %u)\n",
			    __FILE__, __LINE__, __func__,
			    aTr->tres[aEdges[0]].ps->subtreesScore
			    + aTr->tres[aEdges[0]].ps->nodeScore,
			    aTr->tres[aEdges[0]].ps->subtreesScore,
			    + aTr->tres[aEdges[0]].ps->nodeScore,
			    aTr->tres[aEdges[i]].ps->subtreesScore
			    + aTr->tres[aEdges[i]].ps->nodeScore,
			    aTr->tres[aEdges[i]].ps->subtreesScore,
			    aTr->tres[aEdges[i]].ps->nodeScore);
		    abort();
		}
	    }
	}
#endif
    }

    rVal = false;
    RETURN:
    return rVal;
}

// Hold a tree.  If aMaxHeld is exceeded, the tree is not held.  This introduces
// a bias in which trees are held.  There exist algorithms for making this an
// unbiased process, but there is no need for that functionality at the moment.
CxmpInline bool
CxpTrHold(CxtTr *aTr, uint32_t aMaxHold, uint32_t aNeighbor, uint32_t aScore)
{
    bool rVal;

    if (aTr->nheld < aMaxHold)
    {
	CxtTrh *trh;

	// Make sure there is space to store another held tree.
	if (aTr->held == NULL)
	{
	    // Allocate.
	    aTr->held = (CxtTrh *) CxmMalloc(sizeof(CxtTrh));
	    aTr->heldlen = 1;
	}
	else if (aTr->nheld == aTr->heldlen)
	{
	    // Reallocate.
	    aTr->held = (CxtTrh *) CxmRealloc(aTr->held,
					      sizeof(CxtTrh)
					      * aTr->heldlen * 2);
	    aTr->heldlen *= 2;
	}

	// Hold this tree.
	trh = &aTr->held[aTr->nheld];
	trh->neighbor = aNeighbor;
	trh->score = aScore;

	aTr->nheld++;

	rVal = false;
    }
    else
    {
	rVal = true;
    }

    return rVal;
}

// Calculate the Fitch parsimony scores for all TBR neighbors of aTr, and hold
// results according to the function parameters.
CxmpInline void
CxpTrTbrNeighborsMp(CxtTr *aTr, uint32_t aMaxHold, uint32_t aMaxscore,
		    CxtTrHoldHow aHow)
{
    uint32_t neighbor, i, j, k, curMax, score;
    CxtTrEdge bisect, edgeA, edgeB;
    CxtTrNode nodeA, nodeB;
    CxtTrPs *psA, *psB;

    CxmDassert(CxpTrValidate(aTr));

    curMax = aMaxscore;

    // Set up tree holding data structures.
    aTr->nheld = 0;

    // Iteratively (logically) bisect at each edge in the tree.
    neighbor = 0;
    for (i = 0; i < aTr->nedges; i++)
    {
	bisect = aTr->trt[i].bisectEdge;

	// Determine which edges are in each subtree.
	CxpTrBedgesGen(aTr, bisect, &nodeA, &nodeB);

	// Calculate the partial score for each edge in the edge lists.  Don't
	// bother scoring the trees if either subtree exceeds the max score.
	if (CxpTrBisectionEdgeListMp(aTr, aTr->bedges,
				     aTr->nbedgesA, bisect, curMax)
	    || CxpTrBisectionEdgeListMp(aTr, &aTr->bedges[aTr->nbedgesA],
					aTr->nbedgesB, bisect, curMax))
	{
	    neighbor += (aTr->trt[i + 1].offset - aTr->trt[i].offset);
	    continue;
	}

	// Iteratively (logically) reconnect every legitimate pairing of edges
	// between the two subtrees and calculate final parsimony scores.
	for (j = 0; j < aTr->nbedgesA; j++)
	{
	    edgeA = aTr->bedges[j];
	    if (edgeA != CxmTrEdgeNone)
	    {
		psA = aTr->tres[edgeA].ps;
	    }
	    else
	    {
		psA = aTr->trrs[CxmQliFirst(&aTr->trns[nodeA].rings)].ps;
	    }

	    for (k = 0; k < aTr->nbedgesB; k++)
	    {
		// Skip this iteration if the reconnection would result in
		// reversing the bisection.
		if (j == 0 && k == 0)
		{
		    continue;
		}

		edgeB = aTr->bedges[aTr->nbedgesA + k];
		if (edgeB != CxmTrEdgeNone)
		{
		    psB = aTr->tres[edgeB].ps;
		}
		else
		{
		    psB = aTr->trrs[CxmQliFirst(&aTr->trns[nodeB].rings)].ps;
		}

		// Calculate the final parsimony score for this reconnection.
		score = CxpTrMpFscore(aTr, psA, psB, curMax);

		// Hold the tree, if appropriate.
		switch (aHow)
		{
		    case CxeTrHoldBest:
		    {
			if (score == curMax)
			{
			    // This tree is as good as those currently held.
			    CxpTrHold(aTr, aMaxHold, neighbor, score);
			}
			else if (score < curMax)
			{
			    // No trees held, or this tree is better than those
			    // currently held.  Throw away previously held trees
			    // and hold this one.
			    aTr->nheld = 0;
			    CxpTrHold(aTr, aMaxHold, neighbor, score);
			    curMax = score;
			}
			break;
		    }
		    case CxeTrHoldBetter:
		    {
			if (score <= curMax)
			{
			    // No trees held, or this (neighboring) tree is
			    // better than the tree whose neighbors are being
			    // evaluated.
			    CxpTrHold(aTr, aMaxHold, neighbor, score);
			}
			break;
		    }
		    case CxeTrHoldAll:
		    {
			// Hold all trees.
			CxpTrHold(aTr, aMaxHold, neighbor, score);
			break;
		    }
		    default:
		    {
			CxmNotReached();
		    }
		}

		// Increment the neighbor index here.  Due to the possibility of
		// loop short-circuting above, this must happen at the end of
		// the loop body, rather than in the 'for' statement.
		neighbor++;
	    }
	}
    }
}

CxtTr *
CxTrNew(void)
{
    CxtTr *rVal;

    rVal = (CxtTr *) CxmMalloc(sizeof(CxtTr));
    CxpTrNew(rVal);

    return rVal;
}

void
CxTrDelete(CxtTr *aTr)
{
    CxmCheckPtr(aTr);
    CxmAssert(aTr->magic == CxmTrMagic);

    if (aTr->held != NULL)
    {
	CxmFree(aTr->held);
    }

    // This assumes that all nodes are deallocated before CxTrDelete() is
    // called.
    if (aTr->trns != NULL)
    {
	CxmFree(aTr->trns);
    }

    if (aTr->trt != NULL)
    {
	CxmFree(aTr->trt);
    }

    if (aTr->bedges != NULL)
    {
	CxmFree(aTr->bedges);
    }

    // This assumes that all edges are deallocated before CxTrDelete() is
    // called.
    if (aTr->tres != NULL)
    {
	CxmFree(aTr->tres);
	CxmFree(aTr->trrs);
    }

    CxmFree(aTr);
}

uint32_t
CxTrNtaxaGet(CxtTr *aTr)
{
    CxpTrUpdate(aTr);
    CxmDassert(CxpTrValidate(aTr));

    return aTr->ntaxa;
}

uint32_t
CxTrNedgesGet(CxtTr *aTr)
{
    CxpTrUpdate(aTr);
    CxmDassert(CxpTrValidate(aTr));

    return aTr->nedges;
}

CxtTrNode
CxTrBaseGet(CxtTr *aTr)
{
    CxmCheckPtr(aTr);
    CxmAssert(aTr->magic == CxmTrMagic);

    return aTr->base;
}

void
CxTrBaseSet(CxtTr *aTr, CxtTrNode aBase)
{
    CxmCheckPtr(aTr);
    CxmAssert(aTr->magic == CxmTrMagic);
#ifdef CxmDebug
    if (aBase != CxmTrNodeNone)
    {
	CxmDassert(CxpTrNodeValidate(aTr, aBase));
    }
#endif

    aTr->base = aBase;

    aTr->modified = true;
}

void *
CxTrAuxGet(CxtTr *aTr)
{
    CxmCheckPtr(aTr);
    CxmAssert(aTr->magic == CxmTrMagic);

    return aTr->aux;
}

void
CxTrAuxSet(CxtTr *aTr, void *aAux)
{
    CxmCheckPtr(aTr);
    CxmAssert(aTr->magic == CxmTrMagic);

    aTr->aux = aAux;
}

void
CxTrMpFinish(CxtTr *aTr)
{

    CxpTrUpdate(aTr);
    CxmDassert(CxpTrValidate(aTr));

    if (aTr->base != CxmTrNodeNone)
    {
	CxtTrn *trn;
	CxtTrRing ring;
	CxtTre *tre;

	// Clean up the tree.
	trn = &aTr->trns[aTr->base];
	CxmQliForeach(ring, &trn->rings, aTr->trrs, link)
	    {
		// Clean up edge before recursing.
		tre = &aTr->tres[CxTrRingEdgeGet(aTr, ring)];
		if (tre->ps != NULL)
		{
		    CxpTrPsDelete(aTr, tre->ps);
		    tre->ps = NULL;
		}

		// Clean up ring.
		CxpTrMpRingFinish(aTr, ring);

		// Recurse.
		CxpTrMpFinishRecurse(aTr, CxTrRingOtherGet(aTr, ring));
	    }
    }
}

uint32_t
CxTrMpScore(CxtTr *aTr)
{
    uint32_t rVal;
    CxtTrRing ring;

    CxmDassert(CxpTrValidate(aTr));

    if (aTr->base != CxmTrNodeNone
	&& (ring = CxmQliFirst(&aTr->trns[aTr->base].rings)) != CxmTrRingNone)
    {
	CxtTrEdge edge;
	CxtTrPs *psA, *psB;

	edge = CxTrRingEdgeGet(aTr, ring);

	// Calculate partial scores for the subtrees on each end of edge.
	psA = CxpTrMpScoreRecurse(aTr,
				  CxTrEdgeRingGet(aTr, edge, 0),
				  CxmTrEdgeNone);
	psB = CxpTrMpScoreRecurse(aTr,
				  CxTrEdgeRingGet(aTr, edge, 1),
				  CxmTrEdgeNone);

	// Calculate the final score.
	rVal = CxpTrMpFscore(aTr, psA, psB, UINT_MAX);
    }
    else
    {
	rVal = 0;
    }

    return rVal;
}

void
CxTrTbrBestNeighborsMp(CxtTr *aTr, uint32_t aMaxHold)
{
    CxmDassert(CxpTrValidate(aTr));

    CxpTrTbrNeighborsMp(aTr, aMaxHold,
			CxmTrMaxscoreNone,
			CxeTrHoldBest);
}

void
CxTrTbrBetterNeighborsMp(CxtTr *aTr, uint32_t aMaxHold)
{
    uint32_t score;

    CxmDassert(CxpTrValidate(aTr));

    score = CxTrMpScore(aTr);
    CxpTrTbrNeighborsMp(aTr, aMaxHold,
			score > 0 ? score - 1 : 0,
			CxeTrHoldBetter);
}

void
CxTrTbrAllNeighborsMp(CxtTr *aTr)
{
    CxmDassert(CxpTrValidate(aTr));

    CxpTrTbrNeighborsMp(aTr, CxmTrHoldAll,
			CxmTrMaxscoreNone,
			CxeTrHoldAll);
}

void
CxTrHeldFinish(CxtTr *aTr)
{
    CxmDassert(CxpTrValidate(aTr));

    if (aTr->held != NULL)
    {
	CxmFree(aTr->held);
	aTr->held = NULL;
	aTr->heldlen = 0;
	aTr->nheld = 0;
    }
}

uint32_t
CxTrNheldGet(CxtTr *aTr)
{
    CxpTrUpdate(aTr);
    CxmDassert(CxpTrValidate(aTr));

    return aTr->nheld;
}

void
CxTrHeldGet(CxtTr *aTr, uint32_t aHeld, uint32_t *rNeighbor,
	    uint32_t *rScore)
{
    CxtTrh *trh;

    CxpTrUpdate(aTr);
    CxmDassert(CxpTrValidate(aTr));
    CxmCheckPtr(aTr->held);
    CxmAssert(aHeld < aTr->nheld);

    trh = &aTr->held[aHeld];
    *rNeighbor = trh->neighbor;
    *rScore = trh->score;
}
